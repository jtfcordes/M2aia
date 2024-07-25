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

// Qmitk
#include "m2Position.h"

// Qt
#include <QMessageBox>

// mitk image
#include <mitkImage.h>
#include <mitkImageAccessByItk.h>
#include <mitkNodePredicateDataType.h>
#include <mitkPointSet.h>
#include <mitkRotationOperation.h>
#include <mitkScaleOperation.h>
//#include <mitk>
#include <QShortcut>
#include <m2ImzMLSpectrumImage.h>
#include <mitkApplyTransformMatrixOperation.h>
#include <mitkInteractionConst.h>
#include <qpushbutton.h>

const std::string m2Position::VIEW_ID = "org.mitk.views.m2.Position";

void m2Position::SetFocus() {}

void m2Position::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi(parent);

  auto rotateLeft = [this]() { this->Rotate(-0.1); };
  auto rotateRight = [this]() { this->Rotate(0.1); };

  connect(m_Controls.btnPlus1, &QPushButton::clicked, [this]() { this->Rotate(1); });
  connect(m_Controls.btnPlus5, &QPushButton::clicked, [this]() { this->Rotate(5); });
  connect(m_Controls.btnPlus15, &QPushButton::clicked, [this]() { this->Rotate(15); });
  connect(m_Controls.btnPlus45, &QPushButton::clicked, [this]() { this->Rotate(45); });
  connect(m_Controls.btnPlus90, &QPushButton::clicked, [this]() { this->Rotate(90); });
  
  connect(m_Controls.btnMinus1, &QPushButton::clicked, [this]() { this->Rotate(-1); });
  connect(m_Controls.btnMinus5, &QPushButton::clicked, [this]() { this->Rotate(-5); });
  connect(m_Controls.btnMinus15, &QPushButton::clicked, [this]() { this->Rotate(-15); });
  connect(m_Controls.btnMinus45, &QPushButton::clicked, [this]() { this->Rotate(-45); });
  connect(m_Controls.btnMinus90, &QPushButton::clicked, [this]() { this->Rotate(-90); });
  // connect(m_Controls.mirrorH, &QPushButton::clicked, [this]() {this->Mirror(0); });
  // connect(m_Controls.mirrorV, &QPushButton::clicked, [this]() {this->Mirror(1); });

  QShortcut *scRotateLeft = new QShortcut(QKeySequence(Qt::Key_7), parent);
  QShortcut *scRotateRight = new QShortcut(QKeySequence(Qt::Key_9), parent);
  connect(scRotateLeft, &QShortcut::activated, rotateLeft);
  connect(scRotateRight, &QShortcut::activated, rotateRight);

  auto left = [this]() { this->Move({-m_Controls.spnBxStepWidth->value(), 0}); };
  auto right = [this]() { this->Move({m_Controls.spnBxStepWidth->value(), 0}); };
  auto up = [this]() { this->Move({0, -m_Controls.spnBxStepWidth->value()}); };
  auto down = [this]() { this->Move({0, m_Controls.spnBxStepWidth->value()}); };

  connect(m_Controls.btnLeft, &QPushButton::clicked, left);
  connect(m_Controls.btnRight, &QPushButton::clicked, right);
  connect(m_Controls.btnUp, &QPushButton::clicked, up);
  connect(m_Controls.btnDown, &QPushButton::clicked, down);

  QShortcut *scUp = new QShortcut(QKeySequence(Qt::Key_8), parent);
  QShortcut *scDown = new QShortcut(QKeySequence(Qt::Key_2), parent);
  QShortcut *scLeft = new QShortcut(QKeySequence(Qt::Key_4), parent);
  QShortcut *scRight = new QShortcut(QKeySequence(Qt::Key_6), parent);

  connect(scUp, &QShortcut::activated, up);
  connect(scDown, &QShortcut::activated, down);
  connect(scLeft, &QShortcut::activated, left);
  connect(scRight, &QShortcut::activated, right);
}

void m2Position::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*source*/,
                                    const QList<mitk::DataNode::Pointer> &nodes)
{
  // iterate all selected objects, adjust warning visibility
  foreach (mitk::DataNode::Pointer node, nodes)
  {
    if (node.IsNotNull() && dynamic_cast<mitk::Image *>(node->GetData()))
    {
      return;
    }
  }
}

void m2Position::Move(std::array<int, 2> &&vec)
{
  QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
  if (nodes.empty())
    return;

  for (auto node : nodes)
  {
    mitk::BaseData *data = node->GetData();
    if (data)
    {
      // test if this data item is an image or not (could also be a surface or something totally different)
      mitk::Vector3D v;
      v[0] = vec[0] * 10e-4;
      v[1] = vec[1] * 10e-4;
      v[2] = 0;
      if (m2::SpectrumImage *image = dynamic_cast<m2::SpectrumImage *>(data))
      {

        image->ApplyMoveOriginOperation(v);
      }
      else if (mitk::Image *image = dynamic_cast<mitk::Image *>(data))
      {
        auto geometry = image->GetGeometry();
        geometry->Translate(v);
      }
      RequestRenderWindowUpdate();

      // RequestRenderWindowUpdate();
      // auto deriv = GetDataStorage()->GetDerivations(node, mitk::TNodePredicateDataType<mitk::PointSet>::New());
      // if (deriv->Size())
      // {
      //   for (auto p : *deriv)
      //   {
      //     auto geometry = p->GetData()->GetGeometry();
      //     auto pos = geometry->GetOrigin();
      //     auto space = geometry->GetSpacing();
      //     pos[0] = pos[0] + vec.at(0) * space[0];
      //     pos[1] = pos[1] + vec.at(1) * space[1];
      //     geometry->SetOrigin(pos);
      //   }
      // }
    }
  }
}

void m2Position::Mirror(int w)
{
  QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
  if (nodes.empty())
    return;

  for (auto node : nodes)
  {
    mitk::BaseData *data = node->GetData();
    if (data)
    {
      // test if this data item is an image or not (could also be a surface or something totally different)
      m2::SpectrumImage *image = dynamic_cast<m2::SpectrumImage *>(data);
      if (image)
      {
        std::unique_ptr<mitk::ApplyTransformMatrixOperation> op;

        vtkSmartPointer<vtkMatrix4x4> m = vtkMatrix4x4::New();
        m->Identity();
        m->SetElement(w, w, -1.0);
        op.reset(new mitk::ApplyTransformMatrixOperation(
          mitk::EOperations::OpAPPLYTRANSFORMMATRIX, m, image->GetGeometry()->GetCenter()));

        image->ApplyGeometryOperation(op.get());
        RequestRenderWindowUpdate();
      }
    }
  }
}

void m2Position::Rotate(float angleDeg)
{
  QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
  if (nodes.empty())
    return;

  for (auto node : nodes)
  {
    mitk::BaseData *data = node->GetData();
    if (data)
    {
      // test if this data item is an image or not (could also be a surface or something totally different)

      std::unique_ptr<mitk::RotationOperation> op;
      mitk::ScalarType rotAx[3] = {0, 0, 1};

      if (m2::SpectrumImage *image = dynamic_cast<m2::SpectrumImage *>(data))
      {
        op.reset(new mitk::RotationOperation(
          mitk::EOperations::OpROTATE, image->GetGeometry()->GetCenter(), mitk::Vector3D(rotAx), angleDeg));
        image->ApplyGeometryOperation(op.get());
        RequestRenderWindowUpdate();

        // auto deriv = GetDataStorage()->GetDerivations(node, mitk::TNodePredicateDataType<mitk::PointSet>::New());

        // for (auto p : *deriv)
        // {
        //   auto manipulatedGeometry = p->GetData()->GetGeometry()->Clone();
        //   op.reset(new mitk::RotationOperation(mitk::EOperations::OpROTATE,
        //                                        image->GetGeometry()->GetCenter(),
        //                                        mitk::Vector3D(rotAx),
        //                                        angleDeg));
        //   manipulatedGeometry->ExecuteOperation(op.get());
        //   p->GetData()->GetGeometry()->SetIdentity();
        //   p->GetData()->GetGeometry()->Compose(manipulatedGeometry->GetIndexToWorldTransform());
        // }
        // RequestRenderWindowUpdate();
      }
      else if (auto image = dynamic_cast<mitk::Image *>(data))
      {
        op.reset(new mitk::RotationOperation(
          mitk::EOperations::OpROTATE, image->GetGeometry()->GetCenter(), mitk::Vector3D(rotAx), angleDeg));
        auto manipulated = image->GetGeometry()->Clone();
        manipulated->ExecuteOperation(op.get());
        image->GetGeometry()->SetIdentity();
        image->GetGeometry()->Compose(manipulated->GetIndexToWorldTransform());
        RequestRenderWindowUpdate();
      }
    }
  }
}
