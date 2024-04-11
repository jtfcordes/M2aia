/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkSimCLRView.h"

#include <QMessageBox>
#include <m2SpectrumImage.h>
#include <m2SpectrumImageHelper.h>
#include <mitkDockerHelper.h>
#include <mitkNodePredicateFunction.h>
#include <mitkProgressBar.h>
#include <m2IntervalVector.h>

// Don't forget to initialize the VIEW_ID.
const std::string QmitkSimCLRView::VIEW_ID = "org.mitk.views.m2.docker.simclr";

void QmitkSimCLRView::CreateQtPartControl(QWidget *parent)
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

void QmitkSimCLRView::SetFocus()
{
  m_Controls.btnRun->setFocus();
}

void QmitkSimCLRView::EnableWidgets(bool enable)
{
  m_Controls.btnRun->setEnabled(enable);
}

void QmitkSimCLRView::OnStartDockerProcessing()
{
  if (mitk::DockerHelper::CanRunDocker())
  {
    for (auto centroidNode : m_Controls.centroidsSeceltion->GetSelectedNodesStdVector())
      for (auto imageNode : m_Controls.imageSelection->GetSelectedNodesStdVector())
      {

        m2::IntervalVector * centroids = dynamic_cast<m2::IntervalVector *>(centroidNode->GetData());
        if((unsigned int)(m_Controls.batch_size->value()) > centroids->GetIntervals().size()){
          QMessageBox::critical(nullptr, QString("Batch size error") ,QString("Number of centroids must be smaller than the number of samples in a batch."));
          return;
        }

        try
        {
          auto refImage = dynamic_cast<mitk::Image *>(imageNode->GetData());

          mitk::ProgressBar::GetInstance()->AddStepsToDo(3);
          
          mitk::DockerHelper helper("ghcr.io/m2aia/simclr:latest");
          m2::SpectrumImageHelper::AddArguments(helper);

          helper.EnableAutoRemoveContainer(true);
          helper.EnableGPUs(true);

          helper.AddAutoSaveData(imageNode->GetData(), "--imzml", "processData",".imzML");
          helper.AddAutoSaveData(centroidNode->GetData(), "--centroids", "input",".centroids");
          helper.AddApplicationArgument("--model", "/tmp/model.simclr");
          helper.AddApplicationArgument("--epochs", m_Controls.epochs->text().toStdString());
          helper.AddApplicationArgument("--batch_size", m_Controls.batch_size->text().toStdString());
          helper.AddApplicationArgument("--num_clusters", m_Controls.num_clusters->text().toStdString());
          helper.AddApplicationArgument("--num_neighbors", m_Controls.num_neighbors->text().toStdString());
          helper.AddAutoLoadOutput("--umap", "clusters_umap.png");
          helper.AddAutoLoadOutput("--out", "cluster_results.centroids");

          // from data view?
          // helper.AddAutoLoadOutputFolder("--")

          mitk::ProgressBar::GetInstance()->Progress();
          const auto results = helper.GetResults();
          mitk::ProgressBar::GetInstance()->Progress();

          auto resultsIt = std::begin(results);

          auto newNode = mitk::DataNode::New();
          newNode->SetData(*resultsIt);
          newNode->SetName(imageNode->GetName() + "_umap");

          if(auto image = dynamic_cast<mitk::Image *>((*resultsIt).GetPointer())){
            image->GetGeometry()->SetSpacing(refImage->GetGeometry()->GetSpacing());
            image->GetGeometry()->SetOrigin(refImage->GetGeometry()->GetOrigin());
          }

          GetDataStorage()->Add(newNode, const_cast<mitk::DataNode *>(imageNode.GetPointer()));
          ++resultsIt;

          unsigned int i = 0;
          while(resultsIt != std::end(results)){
            auto newNode = mitk::DataNode::New();
            newNode->SetData(*resultsIt);
            newNode->SetName(imageNode->GetName() + "_" + std::to_string(i));
            GetDataStorage()->Add(newNode, const_cast<mitk::DataNode *>(imageNode.GetPointer()));
            m2::CopyNodeProperties(imageNode, newNode);
            ++resultsIt;
            ++i;
          }

          mitk::ProgressBar::GetInstance()->Progress();
        }
        catch (std::exception &e)
        {
          mitk::ProgressBar::GetInstance()->Reset();
          mitkThrow() << "Running SimCLR failed.\nException: " << e.what();
        }
      }
  }
}
