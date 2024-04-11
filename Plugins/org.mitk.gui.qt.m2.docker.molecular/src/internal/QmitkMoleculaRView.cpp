/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkMoleculaRView.h"

#include <QMessageBox>
#include <m2IntervalVector.h>
#include <m2SpectrumImage.h>
#include <mitkCoreServices.h>
#include <mitkDockerHelper.h>
#include <mitkIPreferences.h>
#include <mitkIPreferencesService.h>
#include <mitkImage.h>
#include <mitkNodePredicateFunction.h>
#include <mitkProgressBar.h>
#include <usModuleRegistry.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>

// Don't forget to initialize the VIEW_ID.
const std::string QmitkMoleculaRView::VIEW_ID = "org.mitk.views.m2.docker.molecular";

void QmitkMoleculaRView::CreateQtPartControl(QWidget *parent)
{
  // Setting up the UI is a true pleasure when using .ui files, isn't it?
  m_Controls.setupUi(parent);

  auto NodePredicateIsSpectrumImage = mitk::NodePredicateFunction::New(
    [this](const mitk::DataNode *node) -> bool
    {
      if (auto image = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
        return image->GetImageAccessInitialized();
      return false;
    });

  // m_Controls.imageSelection->SetAutoSelectNewNodes(true);
  m_Controls.imageSelection->SetDataStorage(this->GetDataStorage());
  m_Controls.imageSelection->SetSelectionIsOptional(true);
  m_Controls.imageSelection->SetEmptyInfo(QStringLiteral("Select an image"));
  m_Controls.imageSelection->SetNodePredicate(NodePredicateIsSpectrumImage);

  auto NodePredicateIsCentroidList = mitk::NodePredicateFunction::New(
    [this](const mitk::DataNode *node) -> bool
    {
      if (auto centroidList = dynamic_cast<m2::IntervalVector *>(node->GetData()))
        return ((unsigned int)(centroidList->GetType())) & ((unsigned int)(m2::SpectrumFormat::Centroid));
      return false;
    });
  m_Controls.centroidsSelection->SetDataStorage(this->GetDataStorage());
  m_Controls.centroidsSelection->SetSelectionIsOptional(true);
  m_Controls.centroidsSelection->SetEmptyInfo(QStringLiteral("Select an centroid list"));
  m_Controls.centroidsSelection->SetNodePredicate(NodePredicateIsCentroidList);

  // Wire up the UI widgets with our functionality.
  connect(m_Controls.btnRunMoleculaR, SIGNAL(clicked()), this, SLOT(OnStartMoleculaR()));
}

void QmitkMoleculaRView::SetFocus()
{
  m_Controls.btnRunMoleculaR->setFocus();
}

void QmitkMoleculaRView::EnableWidgets(bool enable)
{
  m_Controls.btnRunMoleculaR->setEnabled(enable);
}

void QmitkMoleculaRView::OnStartMoleculaR()
{
  if (mitk::DockerHelper::CanRunDocker())
  {
    for (auto node : m_Controls.imageSelection->GetSelectedNodes())
    {
      if (auto image = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
      {
        try
        {
          mitk::ProgressBar::GetInstance()->AddStepsToDo(4);
          std::vector<mitk::BaseData::Pointer> data;

          auto *preferencesService = mitk::CoreServices::GetPreferencesService();
          auto *preferences = preferencesService->GetSystemPreferences();
          auto tolerance = preferences->GetFloat("m2aia.signal.Tolerance", 75.0);
          auto centroidNodes = m_Controls.centroidsSelection->GetSelectedNodesStdVector();
          // get current shown image as ion image
          if (centroidNodes.empty())
          {
            auto newImage = mitk::Image::New();
            newImage->Initialize(image);
            mitk::ImagePixelReadAccessor<m2::DisplayImagePixelType, 3> racc(image);
            mitk::ImagePixelWriteAccessor<m2::DisplayImagePixelType, 3> wacc(newImage);
            auto d = image->GetDimensions();
            std::copy(racc.GetData(), racc.GetData() + (d[0]*d[1]*d[2]), wacc.GetData());
            data.push_back(newImage);
          }
          // create ion image for list of centroids
          else
          {
            for (auto centroidNode : centroidNodes)
            {
              if (auto centroids = dynamic_cast<m2::IntervalVector *>(centroidNode->GetData()))
              {
                for (auto i : centroids->GetIntervals())
                {
                  auto newImage = mitk::Image::New();
                  newImage->Initialize(image);
                  image->GetImage(i.x.mean(), tolerance, nullptr, newImage);
                  data.push_back(newImage);
                }
              }
            }
            MITK_INFO << "Size of ionimage vector " << data.size();
          }

          mitk::DockerHelper helper("ghcr.io/m2aia/molecular");
          helper.EnableAutoRemoveContainer(true);
          helper.AddAutoSaveData(data, "--ionimage", "ionimages/image%1%", ".nrrd");
          helper.AddAutoSaveData(image->GetMaskImage(), "--mask", "mask", ".nrrd");
          helper.AddApplicationArgument("--pval", std::to_string(m_Controls.spnBxMPMPValue->value()));
          helper.AddAutoLoadOutput("--out", "mpm.nrrd");
          mitk::ProgressBar::GetInstance()->Progress();

          const auto results = helper.GetResults();
          auto lsImage = mitk::LabelSetImage::New();
          lsImage->InitializeByLabeledImage(dynamic_cast<mitk::Image *>(results[0].GetPointer()));
          mitk::ProgressBar::GetInstance()->Progress();

          // {
          //   using namespace std;
          //   mitk::ImagePixelReadAccessor<mitk::LabelSetImage::PixelType> mAcc(image->GetMaskImage());
          //   mitk::ImagePixelWriteAccessor<mitk::LabelSetImage::PixelType> lsAcc(lsImage);
          //   auto N = accumulate(lsImage->GetDimensions(), lsImage->GetDimensions() + 3, 1, multiplies<unsigned
          //   long>()); transform(mAcc.GetData(), mAcc.GetData() + N, lsAcc.GetData(), lsAcc.GetData(), [](auto a,
          //   auto b) {
          //     return a > 0 ? b : 0;
          //   });
          // }
          mitk::ProgressBar::GetInstance()->Progress();

          auto newNode = mitk::DataNode::New();
          newNode->SetData(lsImage);
          newNode->SetName(node->GetName() + "_mpm_" + std::to_string(image->GetCurrentX()));
          GetDataStorage()->Add(newNode, node);
          lsImage->Update();
          RequestRenderWindowUpdate();

          mitk::Color h;
          h.Set(245 / 255.0, 93 / 255.0, 80 / 255.0);
          lsImage->GetLabel(1)->SetName("Hotspots");
          lsImage->GetLabel(1)->SetColor(h);
          lsImage->GetLabel(1)->SetOpacity(0.1);

          mitk::Color c;
          c.Set(66 / 255.0, 151 / 255.0, 160 / 255.0);
          lsImage->GetLabel(2)->SetName("Coldspots");
          lsImage->GetLabel(2)->SetColor(c);
          lsImage->GetLabel(2)->SetOpacity(0.1);

          lsImage->GetActiveLabelSet()->UpdateLookupTable(1);
          lsImage->GetActiveLabelSet()->UpdateLookupTable(2);
          lsImage->GetActiveLabelSet()->SetActiveLabel(1);
          mitk::ProgressBar::GetInstance()->Progress();
        }
        catch (std::exception &e)
        {
          mitk::ProgressBar::GetInstance()->Reset();
          mitkThrow() << "Running MoleculaR failed.\nException: " << e.what();
        }
      }
    }
  }
}
