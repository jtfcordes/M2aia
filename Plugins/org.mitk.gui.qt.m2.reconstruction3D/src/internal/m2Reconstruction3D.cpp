/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/jtfcordes/m2aia for details.

===================================================================*/

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>
#include <berryPlatform.h>

// Qt
#include <QMessageBox>
#include <QThreadPool>

// mitk
#include <mitkNodePredicateDataType.h>
#include <mitkProgressBar.h>

// m2aia
#include "m2Reconstruction3D.h"
#include <m2DataNodePredicates.h>
#include <m2ElxRegistrationHelper.h>
#include <m2ElxUtil.h>
#include <m2ImzMLSpectrumImage.h>
#include <m2NameDialog.h>

const std::string m2Reconstruction3D::VIEW_ID = "org.mitk.views.m2.reconstruction3D";

void m2Reconstruction3D::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Parent = parent;
  m_Controls.setupUi(parent);

  m_Controls.centroidListSelector->SetDataStorage(GetDataStorage());
  m_Controls.centroidListSelector->SetSelectionIsOptional(true);
  m_Controls.centroidListSelector->SetAutoSelectNewNodes(true);
  m_Controls.centroidListSelector->SetEmptyInfo(QString("Centroid Spectrum"));
  m_Controls.centroidListSelector->SetPopUpTitel(QString("Centroid Spectrum"));
  m_Controls.centroidListSelector->SetNodePredicate(
    mitk::NodePredicateAnd::New(m2::DataNodePredicates::NoActiveHelper,
                                m2::DataNodePredicates::IsOverviewSpectrum,
                                m2::DataNodePredicates::IsCentroidSpectrum));

  connect(m_Controls.btnUpdateList, &QPushButton::clicked, this, &m2Reconstruction3D::OnUpdateList);

  connect(m_Controls.btnStartStacking, SIGNAL(clicked()), this, SLOT(OnStartStacking()));

  {
    m_List1 = m_Controls.listWidget;
    m_List1->setDragEnabled(true);
    m_List1->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_List1->setDragDropMode(QAbstractItemView::NoDragDrop);
    // m_List1->setDefaultDropAction(Qt::MoveAction);
    m_List1->setAlternatingRowColors(true);
    m_List1->setSelectionMode(QAbstractItemView::MultiSelection);
    m_List1->setSelectionBehavior(QAbstractItemView::SelectionBehavior::SelectItems);
    m_List1->setSortingEnabled(true);
  }
  {
    m_List2 = m_Controls.listWidget_2;
    m_List2->setDragEnabled(true);
    m_List2->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_List2->setDragDropMode(QAbstractItemView::NoDragDrop);
    // m_List2->setDefaultDropAction(Qt::MoveAction);
    m_List2->setAlternatingRowColors(true);
    m_List2->setSelectionMode(QAbstractItemView::MultiSelection);
    m_List2->setSelectionBehavior(QAbstractItemView::SelectionBehavior::SelectItems);
    m_List2->setSortingEnabled(true);
  }

  connect(m_Controls.btnRemove,
          &QAbstractButton::clicked,
          this,
          [this]()
          {
            for (auto item : m_List1->selectedItems())
            {
              delete m_List1->takeItem(m_List1->row(item));
            }

            for (auto item : m_List2->selectedItems())
            {
              delete m_List2->takeItem(m_List2->row(item));
            }
          });

  connect(m_Controls.btnSwitch,
          &QAbstractButton::clicked,
          this,
          [this]()
          {
            for (auto item : m_List1->selectedItems())
            {
              m_List2->addItem(m_List1->takeItem(m_List1->row(item)));
            }

            for (auto item : m_List2->selectedItems())
            {
              m_List1->addItem(m_List2->takeItem(m_List2->row(item)));
            }
          });
}

m2Reconstruction3D::DataTuple m2Reconstruction3D::GetImageDataById(unsigned int id, QListWidget *listWidget)
{
  auto referenceIndex = listWidget->item(id)->data(Qt::UserRole).toUInt();
  return m_referenceMap.at(referenceIndex);
}

std::shared_ptr<m2::ElxRegistrationHelper> m2Reconstruction3D::RegistrationStep(
  unsigned int fixedId,
  QListWidget *fixedSource,
  std::shared_ptr<m2::ElxRegistrationHelper> fixedTransformer,
  unsigned int movingId,
  QListWidget *movingSource)
{
  auto fixedData = GetImageDataById(fixedId, fixedSource);
  auto movingData = GetImageDataById(movingId, movingSource);
  auto fixedImage = fixedData.image;

  // check if a transformer exist for fixed image and apply
  if (fixedTransformer && !fixedTransformer->GetTransformation().empty())
    fixedImage = fixedTransformer->WarpImage(fixedImage);

  std::vector<std::string> parameters = GetParameters();

  // start of the registration procedure
  auto elxHelper = std::make_shared<m2::ElxRegistrationHelper>();
  elxHelper->SetImageData(fixedImage, movingData.image);
  elxHelper->SetRegistrationParameters(parameters);
  elxHelper->GetRegistration();
  return elxHelper;
};

std::vector<std::string> m2Reconstruction3D::GetParameters()
{
  auto rigidParameters = m_Controls.textRigid->toPlainText().toStdString();
  m2::ElxUtil::ReplaceParameter(
    rigidParameters, "MaximumNumberOfIterations", std::to_string(m_Controls.RigidMaxIters->value()));

  auto deformableParameters = m_Controls.textDeformable->toPlainText().toStdString();
  m2::ElxUtil::ReplaceParameter(
    deformableParameters, "MaximumNumberOfIterations", std::to_string(m_Controls.DeformableMaxIters->value()));

  std::vector<std::string> parameter{rigidParameters, deformableParameters};
  return parameter;
}

void m2Reconstruction3D::OnStartStacking()
{
  // check input data
  const auto numItems = m_List1->count();
  const bool doMultiModalImageRegistration = m_List1->count() == m_List2->count();

  if (!numItems)
  {
    QMessageBox::warning(m_Parent, "No data available!", "The list for 3D reconstruction must not be empty!");
    return;
  }

  // name of the result stack
  std::array<std::string, 2> stackNames;
  std::array<QListWidget *, 2> lists = {m_List1, m_List2};
  for (int i = 0; i < 1 + doMultiModalImageRegistration; ++i)
  {
    auto d = GetImageDataById(0, lists[i]);

    m2NameDialog dialog(m_Parent);
    if (auto bp = d.image->GetProperty("m2aia.xs.selection.center"))
    {
      auto center = dynamic_cast<mitk::DoubleProperty *>(bp.GetPointer())->GetValueAsString();
      dialog.SetName(std::string("stack") + " " + center);
    }
    if (dialog.exec())
    {
      stackNames[i] = dialog.GetName();
    }
    else
    {
      MITK_WARN << "Cancle stacking.";
      return;
    }
  }
  auto stackSize = m_List1->count();
  double spacingZ = m2::MicroMeterToMilliMeter(m_Controls.spinBoxZSpacing->value());
  // prepare stacks
  auto spectrumImageStack1 = m2::SpectrumImageStack::New(stackSize, spacingZ);
  auto spectrumImageStack2 = m2::SpectrumImageStack::New(stackSize, spacingZ);
  /*
   * Two modalities
   * M2 is optional
   * List size = 6
   * if size of M2 is not equal to size of M1; W2-path is ignored.
   * _M1____M2___________
   * |0|   |0|
   * |1|   |1|
   * |2|   |2| Start index
   * |3|   |3|
   * |4|   |4|
   * |5|   |5|
   * '''''''''''''''
   *
   * M1-M1 order: 2-1 1-0 2-3 3-4 4-5  --> W1
   * M2-W1 order: 2-2 1-1 0-0 3-3 4-4 5-5 --> W2
   */

  // prepare workbench
  mitk::ProgressBar::GetInstance()->AddStepsToDo(numItems - 1);

  const bool UseSubsequentOrdering = !m_Controls.chkBxCoRegistrationToSelected->isChecked();
  const auto currentRow = m_List1->currentRow() < 0 ? numItems / 2 : m_List1->currentRow();

  // Initialize by fixed image of stack 1
  {
    auto M1 = GetImageDataById(currentRow, m_List1);
    auto elxHelper = std::make_shared<m2::ElxRegistrationHelper>();
    elxHelper->SetImageData(M1.image, M1.image);
    elxHelper->SetRegistrationParameters({}); // identity
    spectrumImageStack1->Insert(currentRow, elxHelper);
    mitk::ProgressBar::GetInstance()->Progress();
  }

  // Initialize stack 2
  if (doMultiModalImageRegistration)
  {
    auto elxHelper = RegistrationStep(currentRow, m_List1, nullptr, currentRow, m_List2);
    spectrumImageStack2->Insert(currentRow, elxHelper);
    mitk::ProgressBar::GetInstance()->Progress();
  }

  QThreadPool::globalInstance()->start(
    [&]()
    {
      QThreadPool::globalInstance()->start(
        [&]()
        {
          for (int movingId = currentRow - 1; movingId >= 0; --movingId)
          {
            // stack 1
            int fixedId = UseSubsequentOrdering ? movingId + 1 : currentRow;
            auto fixedTransformer = spectrumImageStack1->GetSliceTransformers().at(fixedId);
            auto elxHelper = RegistrationStep(fixedId, m_List1, fixedTransformer, movingId, m_List1);
            spectrumImageStack1->Insert(movingId, elxHelper);

            // stack 2
            if (doMultiModalImageRegistration)
            {
              elxHelper = RegistrationStep(movingId, m_List1, elxHelper, movingId, m_List2);
              spectrumImageStack2->Insert(movingId, elxHelper);
            }
            mitk::ProgressBar::GetInstance()->Progress();
          }
        });

      QThreadPool::globalInstance()->start(
        [&]()
        {
          for (int movingId = currentRow + 1; movingId < numItems; ++movingId)
          {
            // stack 1
            int fixedId = UseSubsequentOrdering ? movingId - 1 : currentRow;
            auto fixedTransformer = spectrumImageStack1->GetSliceTransformers().at(fixedId);
            auto elxHelper = RegistrationStep(fixedId, m_List1, fixedTransformer, movingId, m_List1);
            spectrumImageStack1->Insert(movingId, elxHelper);

            // stack 2
            if (doMultiModalImageRegistration)
            {
              elxHelper = RegistrationStep(movingId, m_List1, elxHelper, movingId, m_List2);
              spectrumImageStack2->Insert(movingId, elxHelper);
            }
            mitk::ProgressBar::GetInstance()->Progress();
          }
        });

    });

  QThreadPool::globalInstance()->waitForDone();

  spectrumImageStack1->InitializeProcessor();
  spectrumImageStack1->InitializeGeometry();

  auto node = mitk::DataNode::New();
  node->SetData(spectrumImageStack1);
  node->SetName(stackNames[0]);
  GetDataStorage()->Add(node);

  if (doMultiModalImageRegistration)
  {
    spectrumImageStack2->InitializeProcessor();
    spectrumImageStack2->InitializeGeometry();

    node = mitk::DataNode::New();
    node->SetData(spectrumImageStack2);
    node->SetName(stackNames[1] + "(2)");
    GetDataStorage()->Add(node);
  }
}

void m2Reconstruction3D::OnUpdateList()
{
  auto all = this->GetDataStorage()->GetAll();
  m_List1->clear();
  m_List2->clear();
  m_referenceMap.clear();
  unsigned int i = 0;
  // iterate all objects in data storage and create a list widged item
  //  unsigned id = 1;
  for (mitk::DataNode::Pointer node : *all)
  {
    if (node.IsNull())
      continue;
    if (auto data = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
    {
      auto res = this->GetDataStorage()->GetDerivations(node, mitk::TNodePredicateDataType<mitk::PointSet>::New());

      auto *item = new QListWidgetItem();
      item->setText((node->GetName()).c_str());
      item->setData(Qt::UserRole, QVariant::fromValue(i));

      DataTuple tuple;
      tuple.image = data;
      tuple.mask = data->GetMaskImage();
      // tuple.points = data->Get

      m_referenceMap[i] = tuple;

      ++i;
      m_List1->addItem(item);
    }
  }
}
