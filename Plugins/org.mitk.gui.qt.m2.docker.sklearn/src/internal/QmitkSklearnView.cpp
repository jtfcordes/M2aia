/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkSklearnView.h"

#include <QMessageBox>
#include <m2SpectrumImage.h>
#include <m2SpectrumImageHelper.h>
#include <mitkDockerHelper.h>
#include <mitkNodePredicateFunction.h>
#include <mitkProgressBar.h>
#include <m2IntervalVector.h>

// Don't forget to initialize the VIEW_ID.
const std::string QmitkSklearnView::VIEW_ID = "org.mitk.views.m2.docker.sklearn";

void QmitkSklearnView::CreateQtPartControl(QWidget *parent)
{
  // Setting up the UI is a true pleasure when using .ui files, isn't it?
  m_Controls.setupUi(parent);

  auto NodePredicateSpectrumImage = mitk::NodePredicateFunction::New(
    [this](const mitk::DataNode *node) -> bool
    {
      if (auto image = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
      {
        return image->GetImageAccessInitialized();
      }
      return false;
    });

  // m_Controls.imageSelection->SetAutoSelectNewNodes(true);
  m_Controls.imageSelection->SetDataStorage(this->GetDataStorage());
  m_Controls.imageSelection->SetSelectionIsOptional(true);
  m_Controls.imageSelection->SetEmptyInfo(QStringLiteral("Select an image"));
  m_Controls.imageSelection->SetNodePredicate(NodePredicateSpectrumImage);

  auto NodePredicateIsCentroidList = mitk::NodePredicateFunction::New(
    [this](const mitk::DataNode *node) -> bool
    {
      if (auto centroidList = dynamic_cast<m2::IntervalVector *>(node->GetData()))
        return ((unsigned int)(centroidList->GetType())) & ((unsigned int)(m2::SpectrumFormat::Centroid));
      return false;
    });
  m_Controls.centroidsSeceltion->SetDataStorage(this->GetDataStorage());
  m_Controls.centroidsSeceltion->SetSelectionIsOptional(true);
  m_Controls.centroidsSeceltion->SetEmptyInfo(QStringLiteral("Select an centroid list"));
  m_Controls.centroidsSeceltion->SetNodePredicate(NodePredicateIsCentroidList);

  // Wire up the UI widgets with our functionality.
  connect(m_Controls.btnRun, SIGNAL(clicked()), this, SLOT(OnStartDockerProcessing()));
}

void QmitkSklearnView::SetFocus()
{
  m_Controls.btnRun->setFocus();
}

void QmitkSklearnView::EnableWidgets(bool enable)
{
  m_Controls.btnRun->setEnabled(enable);
}

void QmitkSklearnView::OnStartDockerProcessing()
{
  if (mitk::DockerHelper::CanRunDocker())
  {
    for (auto centroidNode : m_Controls.centroidsSeceltion->GetSelectedNodesStdVector())
      for (auto imageNode : m_Controls.imageSelection->GetSelectedNodesStdVector())
      {

        // m2::IntervalVector * centroids = dynamic_cast<m2::IntervalVector *>(centroidNode->GetData());

        try
        {
          auto refImage = dynamic_cast<mitk::Image *>(imageNode->GetData());

          mitk::ProgressBar::GetInstance()->AddStepsToDo(3);
          
          mitk::DockerHelper helper("ghcr.io/m2aia/sklearn:latest");
          m2::SpectrumImageHelper::AddArguments(helper);

          helper.EnableAutoRemoveContainer(true);
          helper.EnableGPUs(false);

          helper.AddAutoSaveData(imageNode->GetData(), "--imzml", "processData",".imzML");
          helper.AddAutoSaveData(centroidNode->GetData(), "--centroids", "input",".centroids");
          helper.AddApplicationArgument("--num_clusters", m_Controls.num_clusters->text().toStdString());
          helper.AddAutoLoadOutput("--out", "result.nrrd");
          // helper.AddLoadLaterOutput("--out_cluster_centers", "kmeans.csv");
          

          // from data view?
          // helper.AddAutoLoadOutputFolder("--")

          mitk::ProgressBar::GetInstance()->Progress();
          const auto results = helper.GetResults();
          mitk::ProgressBar::GetInstance()->Progress();

          auto resultsIt = std::begin(results);

          auto lsImage = mitk::LabelSetImage::New();
          lsImage->InitializeByLabeledImage(dynamic_cast<mitk::Image *>(resultsIt->GetPointer()));
          lsImage->GetGeometry()->SetSpacing(refImage->GetGeometry()->GetSpacing());
          lsImage->GetGeometry()->SetOrigin(refImage->GetGeometry()->GetOrigin());

          auto newNode = mitk::DataNode::New();
          newNode->SetData(lsImage);
          newNode->SetName(imageNode->GetName() + "_k_means");

          
          

          GetDataStorage()->Add(newNode, const_cast<mitk::DataNode *>(imageNode.GetPointer()));

          mitk::ProgressBar::GetInstance()->Progress();
        }
        catch (std::exception &e)
        {
          mitk::ProgressBar::GetInstance()->Reset();
          mitkThrow() << "Running Sklearn failed.\nException: " << e.what();
        }
      }
  }
}
