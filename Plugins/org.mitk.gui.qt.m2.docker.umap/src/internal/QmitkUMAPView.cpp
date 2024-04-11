/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkUMAPView.h"

#include <QMessageBox>
#include <m2IntervalVector.h>
#include <m2SpectrumImage.h>
#include <m2SpectrumImageHelper.h>
#include <mitkDockerHelper.h>
#include <mitkNodePredicateFunction.h>
#include <mitkProgressBar.h>

// Don't forget to initialize the VIEW_ID.
const std::string QmitkUMAPView::VIEW_ID = "org.mitk.views.m2.docker.umap";

void QmitkUMAPView::CreateQtPartControl(QWidget *parent)
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
  m_Controls.centroidsSelection->SetDataStorage(this->GetDataStorage());
  m_Controls.centroidsSelection->SetSelectionIsOptional(true);
  m_Controls.centroidsSelection->SetEmptyInfo(QStringLiteral("Select an centroid list"));
  m_Controls.centroidsSelection->SetNodePredicate(NodePredicateIsCentroidList);

  // Wire up the UI widgets with our functionality.
  connect(m_Controls.btnRun, SIGNAL(clicked()), this, SLOT(OnStartDockerProcessing()));
}

void QmitkUMAPView::SetFocus()
{
  m_Controls.btnRun->setFocus();
}

void QmitkUMAPView::EnableWidgets(bool enable)
{
  m_Controls.btnRun->setEnabled(enable);
}

void QmitkUMAPView::OnStartDockerProcessing()
{
  if (mitk::DockerHelper::CanRunDocker())
  {
    for (auto centroidNode : m_Controls.centroidsSelection->GetSelectedNodesStdVector())
      for (auto imageNode : m_Controls.imageSelection->GetSelectedNodesStdVector())
      {
        try
        {

          mitk::ProgressBar::GetInstance()->AddStepsToDo(1);

          // docker image name
          mitk::DockerHelper helper("ghcr.io/m2aia/umap:latest");
          m2::SpectrumImageHelper::AddArguments(helper);

          helper.EnableAutoRemoveContainer(true);
          helper.EnableGPUs(false);

          helper.AddAutoSaveData(imageNode->GetData(), "--imzml", "processData",".imzML");
          helper.AddAutoSaveData(centroidNode->GetData(), "--centroids", "input",".centroids");
          helper.AddAutoLoadOutput("--out", "umap.nrrd");
          helper.AddApplicationArgument("--n_neighbors", m_Controls.n_neighbors->text().toStdString());
          helper.AddApplicationArgument("--n_components", m_Controls.n_components->text().toStdString());
          helper.AddApplicationArgument("--metric", m_Controls.metric->currentText().toStdString());
          helper.AddApplicationArgument("--n_epochs", m_Controls.n_epochs->text().toStdString());
          helper.AddApplicationArgument("--learning_rate", m_Controls.learning_rate->text().toStdString());
          helper.AddApplicationArgument("--min_dist", m_Controls.min_dist->text().toStdString());
          helper.AddApplicationArgument("--spread", m_Controls.spread->text().toStdString());
          helper.AddApplicationArgument("--local_connectivity", m_Controls.local_connectivity->text().toStdString());

          // start processing
          const auto results = helper.GetResults();

          // feed back results into MTIK
          auto refImage = dynamic_cast<mitk::Image *>(imageNode->GetData());
          auto image = dynamic_cast<mitk::Image *>(results[0].GetPointer());
          image->GetGeometry()->SetSpacing(refImage->GetGeometry()->GetSpacing());
          image->GetGeometry()->SetOrigin(refImage->GetGeometry()->GetOrigin());
          auto newNode = mitk::DataNode::New();
          newNode->SetData(image);
          newNode->SetName(imageNode->GetName() + "_umap");

          GetDataStorage()->Add(newNode, const_cast<mitk::DataNode *>(imageNode.GetPointer()));

          mitk::ProgressBar::GetInstance()->Progress();
        }
        catch (std::exception &e)
        {
          mitk::ProgressBar::GetInstance()->Reset();
          mitkThrow() << "Running UMAP failed.\nException: " << e.what();
        }
      }
  }
}
