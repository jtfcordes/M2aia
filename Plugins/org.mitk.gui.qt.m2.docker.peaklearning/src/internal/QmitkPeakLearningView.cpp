/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkPeakLearningView.h"

#include <QMessageBox>
#include <m2SpectrumImage.h>
#include <m2SpectrumImageHelper.h>
#include <mitkDockerHelper.h>
#include <mitkLabelSetImage.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateFunction.h>
#include <mitkProgressBar.h>

// Don't forget to initialize the VIEW_ID.
const std::string QmitkPeakLearningView::VIEW_ID = "org.mitk.views.m2.docker.peaklearning";

void QmitkPeakLearningView::CreateQtPartControl(QWidget *parent)
{
  // Setting up the UI is a true pleasure when using .ui files, isn't it?
  m_Controls.setupUi(parent);

  auto NodePredicateIsContinuousSpectrumImage = mitk::NodePredicateFunction::New(
    [this](const mitk::DataNode *node) -> bool
    {
      if (auto spectrumImage = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
        return ((unsigned int)(spectrumImage->GetSpectrumType().Format)) &
               ((unsigned int)(m2::SpectrumFormat::Continuous));
      return false;
    });

  // m_Controls.imageSelection->SetAutoSelectNewNodes(true);
  m_Controls.imageSelection->SetDataStorage(this->GetDataStorage());
  m_Controls.imageSelection->SetSelectionIsOptional(true);
  m_Controls.imageSelection->SetEmptyInfo(QStringLiteral("Select an image"));
  m_Controls.imageSelection->SetNodePredicate(NodePredicateIsContinuousSpectrumImage);

  m_Controls.maskSelection->SetDataStorage(this->GetDataStorage());
  m_Controls.maskSelection->SetSelectionIsOptional(true);
  m_Controls.maskSelection->SetEmptyInfo(QStringLiteral("Select a mask image"));
  m_Controls.maskSelection->SetNodePredicate(mitk::TNodePredicateDataType<mitk::LabelSetImage>::New());

  // Wire up the UI widgets with our functionality.
  connect(m_Controls.btnRun, SIGNAL(clicked()), this, SLOT(OnStartDockerProcessing()));
}

void QmitkPeakLearningView::SetFocus()
{
  m_Controls.btnRun->setFocus();
}

void QmitkPeakLearningView::EnableWidgets(bool enable)
{
  m_Controls.btnRun->setEnabled(enable);
}

void QmitkPeakLearningView::OnStartDockerProcessing()
{
  if (mitk::DockerHelper::CanRunDocker())
  {
    for (auto node : m_Controls.imageSelection->GetSelectedNodes())
    {
      if (auto image = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
      {
        if (!image->GetImageAccessInitialized())
          return;

        mitk::Image::Pointer mask;
        if (auto maskNode = m_Controls.maskSelection->GetSelectedNode())
          mask = dynamic_cast<mitk::Image *>(maskNode->GetData());
        else
          mask = image->GetMaskImage();
        
        try
        {
          mitk::ProgressBar::GetInstance()->AddStepsToDo(3);

          mitk::DockerHelper helper("ghcr.io/m2aia/msipl:latest");
          m2::SpectrumImageHelper::AddArguments(helper);
          
          helper.EnableGPUs(true);
          helper.EnableAutoRemoveContainer(true);

          helper.AddAutoSaveData(image, "--imzml", "processData", ".imzML");
          helper.AddAutoSaveData(mask, "--mask", "mask", ".nrrd");
          helper.AddApplicationArgument("--latent_dim", m_Controls.latent_dim->text().toStdString());
          helper.AddApplicationArgument("--interim_dim", m_Controls.interim_dim->text().toStdString());
          helper.AddApplicationArgument("--beta", m_Controls.beta->text().toStdString());
          helper.AddApplicationArgument("--epochs", m_Controls.epochs->text().toStdString());
          helper.AddApplicationArgument("--batch_size", m_Controls.batch_size->text().toStdString());
          helper.AddAutoLoadOutput("--csv", "peaks.csv");
          helper.AddLoadLaterOutput("--model", "model.h5");

          mitk::ProgressBar::GetInstance()->Progress();
          const auto results = helper.GetResults();
          mitk::ProgressBar::GetInstance()->Progress();

          auto newNode = mitk::DataNode::New();
          newNode->SetData(results[0]);
          newNode->SetName(node->GetName() + "_msiPL");
          GetDataStorage()->Add(newNode, node);


          m2::CopyNodeProperties(node, newNode);

          mitk::ProgressBar::GetInstance()->Progress();
        }
        catch (std::exception &e)
        {
          mitk::ProgressBar::GetInstance()->Reset();
          mitkThrow() << "Running PeakLearning failed.\nException: " << e.what();
        }
      }
    }
  }
}
