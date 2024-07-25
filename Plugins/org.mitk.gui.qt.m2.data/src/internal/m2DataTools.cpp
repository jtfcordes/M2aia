/*===================================================================

Mass Spectrometry Imaging applications for interactive
analysis in MITK (M2aia)

Copyright (c) Jonas Cordes, Hochschule Mannheim.
Division of Medical Informatics.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/

#include "m2DataTools.h"

#include <QmitkRenderWindow.h>
#include <m2ImzMLSpectrumImage.h>
#include <m2SpectrumImage.h>
#include <m2UIUtils.h>
#include <mitkLayoutAnnotationRenderer.h>
#include <mitkLookupTableProperty.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>

const std::string m2DataTools::VIEW_ID = "org.mitk.views.m2.DataTools";

void m2DataTools::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi(parent);
  m_Parent = parent;

  {
    m_Controls.ReferenceLevelWindowSelection->SetDataStorage(GetDataStorage());
    m_Controls.ReferenceLevelWindowSelection->SetNodePredicate(
      mitk::NodePredicateAnd::New(mitk::TNodePredicateDataType<mitk::Image>::New(),
                                  mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object"))));
    m_Controls.ReferenceLevelWindowSelection->SetSelectionIsOptional(true);
    m_Controls.ReferenceLevelWindowSelection->SetEmptyInfo(QString("Reference image selection"));
    m_Controls.ReferenceLevelWindowSelection->SetPopUpTitel(QString("Image"));

    m_Controls.ReferenceSelectionForScaleBar->SetDataStorage(GetDataStorage());
    m_Controls.ReferenceSelectionForScaleBar->SetNodePredicate(
      mitk::NodePredicateAnd::New(mitk::TNodePredicateDataType<mitk::Image>::New(),
                                  mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object"))));
    m_Controls.ReferenceSelectionForScaleBar->SetSelectionIsOptional(true);
    m_Controls.ReferenceSelectionForScaleBar->SetEmptyInfo(QString("Reference image selection"));
    m_Controls.ReferenceSelectionForScaleBar->SetPopUpTitel(QString("Image"));
  }

  // disable reference point set
  m_Controls.refPointSetGroup->setVisible(false);

  connect(m_Controls.btnEqualizeLW, &QAbstractButton::clicked, this, &m2DataTools::OnEqualizeLW);
  connect(m_Controls.resetTiling, &QAbstractButton::clicked, this, &m2DataTools::OnResetTiling);
  connect(m_Controls.applyTiling, &QAbstractButton::clicked, this, &m2DataTools::OnApplyTiling);
  m_Controls.ReferenceSelectionForScaleBar->setEnabled(true);

  connect(m_Controls.ScaleBar,
          &QGroupBox::toggled,
          this,
          [this](bool state)
          {
            if (state && !m_Controls.ReferenceSelectionForScaleBar->GetSelectedNode())
            {
              m_Controls.ScaleBar->setChecked(false);
              return;
            }
            // m_Controls.sclaeBarLabel->setEnabled(true);
            m_Controls.ReferenceSelectionForScaleBar->setEnabled(true);
            m_ColorBarAnnotations[0]->SetVisibility(state);
            UpdateColorBarAndRenderWindows();
            RequestRenderWindowUpdate();
          });

  this->m_ColorBarAnnotations.clear();
  for (int i = 0; i < 2; ++i)
  {
    auto cbAnnotation = mitk::ColorBarAnnotation::New();
    this->m_ColorBarAnnotations.push_back(cbAnnotation);
    cbAnnotation->SetFontSize(20);
    cbAnnotation->SetOrientation(1);
    cbAnnotation->SetVisibility(0);

    m_Controls.scaleBarFontSize->setValue(cbAnnotation->GetFontSize());
    m_Controls.scaleBarOrientation->setCurrentIndex(cbAnnotation->GetOrientation());
    // m_Controls.scaleBarLenght->setValue(cbAnnotation->GetLenght());
    // m_Controls.scaleBarWidth->setValue(cbAnnotation->GetWidth());

    connect(m_Controls.scaleBarFontSize,
            &QSlider::sliderMoved,
            this,
            [this, i](int pos)
            {
              m_ColorBarAnnotations[i]->SetFontSize(pos);
              UpdateColorBarAndRenderWindows();
              RequestRenderWindowUpdate();
            });

    connect(m_Controls.scaleBarLenght,
            &QSlider::sliderMoved,
            this,
            [this, i](int pos)
            {
              Q_UNUSED(pos)
              // m_ColorBarAnnotations[i]->SetLength(pos);
              UpdateColorBarAndRenderWindows();
              RequestRenderWindowUpdate();
            });

    connect(m_Controls.scaleBarWidth,
            &QSlider::sliderMoved,
            this,
            [this, i](int pos)
            {
              Q_UNUSED(pos)
              // m_ColorBarAnnotations[i]->SetWidth(pos);
              UpdateColorBarAndRenderWindows();
              RequestRenderWindowUpdate();
            });

    connect(m_Controls.scaleBarOrientation,
            qOverload<int>(&QComboBox::currentIndexChanged),
            this,
            [this, i](int pos)
            {
              m_ColorBarAnnotations[i]->SetOrientation(pos);
              UpdateColorBarAndRenderWindows();
              RequestRenderWindowUpdate();
            });
  }
}

void m2DataTools::OnResetTiling()
{
  auto allNodes = m2::UIUtils::AllNodes(GetDataStorage());

  //	unsigned int maxWidth = 0, maxHeight = 0;
  if (allNodes->Size() == 0)
    return;

  for (auto &e : *allNodes)
  {
    double initP[] = {0, 0, 0};
    mitk::Point3D origin(initP);
    mitk::Point3D prevOrigin(initP);
    if (auto *image = dynamic_cast<m2::SpectrumImage *>(e->GetData()))
    {
      prevOrigin = image->GetGeometry()->GetOrigin();
      origin = image->GetGeometry()->GetOrigin();

      origin[0] = image->GetPropertyValue<double>("[IMS:1000053] absolute position offset x", 0);
      origin[1] = image->GetPropertyValue<double>("[IMS:1000054] absolute position offset y", 0);
      origin[2] = image->GetPropertyValue<double>("absolute position offset z", 0);
      image->GetGeometry()->SetOrigin(origin);
      std::vector<mitk::BaseData *> imageList{image->GetIndexImage(), image->GetMaskImage(), image->GetPoints()};

      for (auto kv : image->GetNormalizationImages())
        imageList.push_back(kv.second.image);
      

      for (auto current : imageList)
        if (current && current->GetGeometry())

          current->GetGeometry()->SetOrigin(origin);
    }

    double dx = origin[0] - prevOrigin[0];
    double dy = origin[1] - prevOrigin[1];

    auto der = this->GetDataStorage()->GetDerivations(e, mitk::TNodePredicateDataType<mitk::PointSet>::New());
    for (auto &&e : *der)
    {
      auto pts = dynamic_cast<mitk::PointSet *>(e->GetData());
      for (auto p = pts->Begin(); p != pts->End(); ++p)
      {
        auto &pp = p->Value();
        pp[0] += dx;
        pp[1] += dy;
      }
    }
  }
  mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
}

void m2DataTools::OnEqualizeLW()
{
  
  auto allNodes = GetDataStorage()->GetSubset(m_Controls.ReferenceLevelWindowSelection->GetNodePredicate());
  // m2::UIUtils::AllNodes(GetDataStorage());

  if (auto node = this->m_Controls.ReferenceLevelWindowSelection->GetSelectedNode())
  {
    mitk::LevelWindow lw_ref;
    node->GetLevelWindow(lw_ref);

    for (auto &n : *allNodes)
      n->SetLevelWindow(lw_ref);
  }

  RequestRenderWindowUpdate();
}

void m2DataTools::OnApplyTiling()
{
  auto rows = m_Controls.mosaicRows->value();
  unsigned int maxWidth = 0, maxHeight = 0;

  // Spectrum Image Base nodes should never be child nodes!
  std::vector<mitk::DataNode::Pointer> nodes;
  auto allNode = GetDataStorage()->GetAll();
  for (auto node : *allNode)
    if (auto image = dynamic_cast<mitk::Image *>(node->GetData()))
    {
      if(!node->IsVisible(nullptr))
        continue;

      maxWidth = std::max(maxWidth, image->GetDimensions()[0]);
      maxHeight = std::max(maxHeight, image->GetDimensions()[1]);
      nodes.push_back(node);
    }

  int nodesInRow = std::ceil(nodes.size() / double(rows));
  MITK_INFO << "NodesInRow: " << nodesInRow;
  if (nodesInRow < 1)
    return;

  std::sort(nodes.begin(),
            nodes.end(),
            [](mitk::DataNode::Pointer &a, mitk::DataNode::Pointer &b) -> bool
            { return a->GetName().compare(b->GetName()) < 0; });

  int i = 0;
  for (auto node : nodes)
  {
    if(!node->IsVisible(nullptr))
      continue;
    // SpectrumImage Nodes
    mitk::Point3D origin, prevOrigin;
    if (auto image = dynamic_cast<mitk::Image *>(node->GetData()))
    {
      prevOrigin = image->GetGeometry()->GetOrigin();
      origin[0] = maxWidth * int(i % nodesInRow) * image->GetGeometry()->GetSpacing()[0];
      origin[1] = maxHeight * int(i / nodesInRow) * image->GetGeometry()->GetSpacing()[1];
      origin[2] = double(0.0);
      image->GetGeometry()->SetOrigin(origin);
      
      if(auto spectrumImage = dynamic_cast<m2::SpectrumImage *>(node->GetData())){
        std::vector<mitk::BaseData *> imageList{spectrumImage->GetIndexImage(), spectrumImage->GetMaskImage(), spectrumImage->GetPoints()};
        for (auto kv : spectrumImage->GetNormalizationImages())
          imageList.push_back(kv.second.image);
      
        for (auto current : imageList)
          if (current && current->GetGeometry())
            current->GetGeometry()->SetOrigin(origin);
      }
      
     
    }
    else
    {
      continue;
    }

    // child node images
    auto predicateImage = mitk::TNodePredicateDataType<mitk::Image>::New();
    auto childNodes = this->GetDataStorage()->GetDerivations(node, predicateImage);
    for (auto child : *childNodes)
      if (auto image = dynamic_cast<mitk::Image *>(child->GetData()))
        image->GetGeometry()->SetOrigin(origin);

    // child node pointsets
    auto predicatePointSet = mitk::TNodePredicateDataType<mitk::PointSet>::New();
    double dx = origin[0] - prevOrigin[0];
    double dy = origin[1] - prevOrigin[1];
    childNodes = this->GetDataStorage()->GetDerivations(node, predicatePointSet);
    for (auto child : *childNodes)
    {
      if (auto pts = dynamic_cast<mitk::PointSet *>(child->GetData()))
      {
        for (auto p = pts->Begin(); p != pts->End(); ++p)
        {
          auto &pp = p->Value();
          pp[0] += dx;
          pp[1] += dy;
        }
      }
    }

    ++i;
  }
  mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
  // this->RequestRenderWindowUpdate();
}

void m2DataTools::UpdateColorBarAndRenderWindows()
{
  mitk::ColorBarAnnotation::Pointer cbAnnotation;
  auto lookuptabel = mitk::LookupTableProperty::New();

  cbAnnotation = m_ColorBarAnnotations[0];

  auto renderer = GetRenderWindowPart()->GetQmitkRenderWindow("axial")->GetRenderer();
  mitk::LayoutAnnotationRenderer::AddAnnotation(cbAnnotation, renderer);

  auto node = this->m_Controls.ReferenceSelectionForScaleBar->GetSelectedNode();
  if (node)
  {
    if (node->GetProperty(lookuptabel, "LookupTable"))
      cbAnnotation->SetLookupTable(lookuptabel->GetValue()->GetVtkLookupTable());
  }
}

void m2DataTools::UpdateLevelWindow(const mitk::DataNode *node)
{
  if (auto msImageBase = dynamic_cast<mitk::Image *>(node->GetData()))
  {
    mitk::LevelWindow lw;
    node->GetLevelWindow(lw);
    lw.SetAuto(msImageBase);
    const_cast<mitk::DataNode *>(node)->SetLevelWindow(lw);
  }
}
