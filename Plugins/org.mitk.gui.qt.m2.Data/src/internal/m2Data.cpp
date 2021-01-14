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

#include "m2Data.h"

#include <QFutureWatcher>
#include <QShortcut>
#include <QmitkMultiNodeSelectionWidget.h>
#include <QmitkRenderWindow.h>
#include <QmitkSingleNodeSelectionWidget.h>
#include <QtConcurrent>
#include <berryPlatformUI.h>
#include <boost/format.hpp>
#include <itkRescaleIntensityImageFilter.h>
#include <itksys/SystemTools.hxx>
#include <mitkCameraController.h>
#include <mitkIOUtil.h>
#include <mitkImage.h>
#include <mitkImageAccessByItk.h>
#include <mitkImageCast.h>
#include <mitkLabelSetImageConverter.h>
#include <mitkLayoutAnnotationRenderer.h>
#include <mitkLookupTableProperty.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkTimer.h>

const std::string m2Data::VIEW_ID = "org.mitk.views.m2.Data";

// 20201023: custom selection service did not work as expected
// void m2Data::SetSelectionProvider()
//{
//  this->GetSite()->SetSelectionProvider(m_SelectionProvider);
//}

void m2Data::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi(parent);

  // 20201023: custom selection service did not work as expected
  // m_CustomSelectionChangedListener.reset(
  //  new berry::SelectionChangedAdapter<m2Data>(this, &m2Data::OnCustomSelectionChanged));
  // GetSite()->GetWorkbenchWindow()->GetSelectionService()->AddSelectionListener(m_CustomSelectionChangedListener.data());

  // 20201023: use communciation service
  auto serviceRef = m2::CommunicationService::Instance();
  connect(serviceRef, SIGNAL(GrabIonImage(qreal, qreal)), this, SLOT(OnGrabIonImage(qreal, qreal)));
  connect(serviceRef,
          SIGNAL(RequestProcessingNodes(const QString &)),
          this,
          SLOT(OnProcessingNodesRequested(const QString &)));

  // 20201023: custom selection service did not work as expected
  //  this->m_SelectionProvider = new m2::SelectionProvider();

  m_Controls.grpBoxCommon->setChecked(false);

  {
    auto m_MassSpecPredicate = mitk::TNodePredicateDataType<m2::MSImageBase>::New();
    m_MassSpecDataNodeSelectionWidget = new QmitkMultiNodeSelectionWidget();
    m_MassSpecDataNodeSelectionWidget->SetDataStorage(GetDataStorage());
    m_MassSpecDataNodeSelectionWidget->SetNodePredicate(
      mitk::NodePredicateAnd::New(mitk::TNodePredicateDataType<m2::MSImageBase>::New(),
                                  mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object"))));
    m_MassSpecDataNodeSelectionWidget->SetSelectionIsOptional(true);
    m_MassSpecDataNodeSelectionWidget->SetEmptyInfo(QString("MS image selection"));
    m_MassSpecDataNodeSelectionWidget->SetPopUpTitel(QString("MS image"));
    ((QVBoxLayout *)(parent->layout()))->insertWidget(1, m_MassSpecDataNodeSelectionWidget);
  }

  connect(m_Controls.btnOpenPointSetWidget, &QAbstractButton::clicked, []() {
    if (auto platform = berry::PlatformUI::GetWorkbench())
    {
      if (auto workbench = platform->GetActiveWorkbenchWindow())
      {
        if (auto page = workbench->GetActivePage())
        {
          if (page.IsNotNull())
          {
            try
            {
              page->ShowView("org.mitk.views.pointsetinteraction", "", 1);
            }
            catch (berry::PartInitException e)
            {
              BERRY_ERROR << "Error: " << e.what() << std::endl;
            }
          }
        }
      }
    }
  });

  connect(m_Controls.btnAddReferencepointSet, &QAbstractButton::clicked, [this]() {
    auto nodes = this->AllNodes();
    for (auto node : *nodes)
    {
      bool pointSetExist = false;
      auto childNodes = this->GetDataStorage()->GetDerivations(node);
      for (auto child : *childNodes)
      {
        if (dynamic_cast<mitk::PointSet *>(child->GetData()))
        {
          pointSetExist = true;
          break;
        }
      }
      if (!pointSetExist)
      {
        auto img = dynamic_cast<m2::MSImageBase *>(node->GetData());
        auto pointSet = mitk::PointSet::New();
        pointSet->SetGeometry(img->GetGeometry());
        auto dataNode = mitk::DataNode::New();
        dataNode->SetData(pointSet);
        dataNode->SetName(node->GetName());
        this->GetDataStorage()->Add(dataNode, node);

        dataNode->SetFloatProperty("pointsize", img->GetGeometry()->GetSpacing()[0] * 3);
        dataNode->SetFloatProperty("point 2D size", img->GetGeometry()->GetSpacing()[0] * 3);
        img->GetImageArtifacts()["references"] = pointSet;
      }
    }
  });

  /* connect(m_Controls.btnUpdate, &QAbstractButton::clicked, [this] {
     auto nodes = this->AllNodes();
     for (auto node : *nodes)
     {
       if (m2::ImzMLMassSpecImage::Pointer imzml = dynamic_cast<m2::ImzMLMassSpecImage *>(node->GetData()))
       {
         auto childs = this->GetDataStorage()->GetDerivations(node, nullptr, true);
         for (auto child : *childs)
           this->GetDataStorage()->Remove(child);
         this->GetDataStorage()->Remove(node);

         imzml->GetImageArtifacts().clear();
         imzml->GetSpectraArtifacts().clear();
         imzml->SetImageGeometryInitialized(false);
         imzml->InitializeGeometry();

         imzml->SetImageAccessInitialized(false);

         auto newNode = mitk::DataNode::New();
         newNode->SetData(imzml);
         newNode->SetName(node->GetName());
         this->GetDataStorage()->Add(newNode);
       }
     }
   });*/

  connect(m_Controls.applyTiling, &QPushButton::clicked, this, &m2Data::OnApplyTiling);
  connect(m_Controls.resetTiling, &QPushButton::clicked, this, &m2Data::OnResetTiling);

  // connect(m_Controls.btnCopyRescale, &QPushButton::clicked, this, &m2Data::OnCopyRescale);

  /*connect(m_Controls.radioButtonMeanSpectrum, &QRadioButton::toggle, this, [this](auto t) {
    if (t) this->m_CurrentOverviewSpectrumType = m2::MSImageBase::OverviewSpectrumType::Mean;
  });*/

  connect(m_Controls.rbtnTolPPM, &QAbstractButton::clicked, this, [&] { m_Controls.spnBxTol->setSuffix(" ppm"); });

  connect(m_Controls.rbtnTolDa, &QAbstractButton::clicked, this, [&] { m_Controls.spnBxTol->setSuffix(" Da"); });

  connect(m_Controls.btnGrabIonImage, &QAbstractButton::clicked, this, [&] {
    OnGrabIonImage(m_Controls.spnBxMz->value(), -1);
  });

  connect(m_Controls.scaleBarEnableIon, &QRadioButton::toggled, this, [this](bool state) {
    m_ColorBarAnnotations[0]->SetVisibility(state);
    UpdateColorBarAndRenderWindows();
    RequestRenderWindowUpdate();
  });

  connect(m_Controls.scaleBarEnableNormalization, &QRadioButton::toggled, this, [this](bool state) {
    m_ColorBarAnnotations[1]->SetVisibility(state);
    UpdateColorBarAndRenderWindows();
    RequestRenderWindowUpdate();
  });

  connect(m_Controls.btnAddIonImageReference, &QPushButton::clicked, this, &m2Data::EmitIonImageReference);

  connect(m_Controls.btnEqualizeLW, &QPushButton::clicked, this, &m2Data::OnEqualizeLW);

  // step through
  QShortcut *shortcutLeft = new QShortcut(QKeySequence(Qt::Key_Left), parent);
  connect(shortcutLeft, SIGNAL(activated()), this, SLOT(OnPrevIonImage()));

  QShortcut *shortcutRight = new QShortcut(QKeySequence(Qt::Key_Right), parent);
  connect(shortcutRight, SIGNAL(activated()), this, SLOT(OnNextIonImage()));

  // add ion image
  QShortcut *shortcutEnter = new QShortcut(QKeySequence(Qt::Key_Return), parent);
  connect(shortcutEnter, SIGNAL(activated()), this, SLOT(EmitIonImageReference()));

  // connect(m_Controls.btnPrev, &QPushButton::clicked, this, &m2Data::OnPrevIonImage);
  // connect(m_Controls.btnNext, &QPushButton::clicked, this, &m2Data::OnNextIonImage);

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
    m_Controls.scaleBarLenght->setValue(cbAnnotation->GetLenght());
    m_Controls.scaleBarWidth->setValue(cbAnnotation->GetWidth());

    connect(m_Controls.scaleBarForceUpdate, &QPushButton::clicked, this, [this] {
      UpdateColorBarAndRenderWindows();
      RequestRenderWindowUpdate();
    });

    connect(m_Controls.scaleBarFontSize, &QSlider::sliderMoved, this, [this, i](int pos) {
      m_ColorBarAnnotations[i]->SetFontSize(pos);
      UpdateColorBarAndRenderWindows();
      RequestRenderWindowUpdate();
    });

    connect(m_Controls.scaleBarLenght, &QSlider::sliderMoved, this, [this, i](int pos) {
      m_ColorBarAnnotations[i]->SetLength(pos);
      UpdateColorBarAndRenderWindows();
      RequestRenderWindowUpdate();
    });

    connect(m_Controls.scaleBarWidth, &QSlider::sliderMoved, this, [this, i](int pos) {
      m_ColorBarAnnotations[i]->SetWidth(pos);
      UpdateColorBarAndRenderWindows();
      RequestRenderWindowUpdate();
    });

    connect(m_Controls.scaleBarOrientation, qOverload<int>(&QComboBox::currentIndexChanged), this, [this, i](int pos) {
      m_ColorBarAnnotations[i]->SetOrientation(pos);
      UpdateColorBarAndRenderWindows();
      RequestRenderWindowUpdate();
    });
  }
}

m2::NormalizationStrategyType m2Data::GuiToNormalizationStrategyType()
{
  // if (this->Controls()->radioButtonMedianNormalization->isChecked())
  //  return m2::NormalizationStrategyType::Median;
  if (this->Controls()->radioButtonTICNormalization->isChecked())
    return m2::NormalizationStrategyType::TIC;
  if (this->Controls()->radioButtonNoneNormalization->isChecked())
    return m2::NormalizationStrategyType::None;
  if (this->Controls()->radioButtonInFileNormalization->isChecked())
    return m2::NormalizationStrategyType::InFile;
  if (this->Controls()->radioButtonSumNormalization->isChecked())
    return m2::NormalizationStrategyType::Sum;
  return m2::NormalizationStrategyType::None;
}

m2::IonImageGrabStrategyType m2Data::GuiToIonImageGrabStrategyType()
{
  if (this->Controls()->radioButtonMax->isChecked())
    return m2::IonImageGrabStrategyType::Maximum;
  if (this->Controls()->radioButtonMean->isChecked())
    return m2::IonImageGrabStrategyType::Mean;
  if (this->Controls()->radioButtonMedian->isChecked())
    return m2::IonImageGrabStrategyType::Median;
  if (this->Controls()->radioButtonSum->isChecked())
    return m2::IonImageGrabStrategyType::Sum;
  return m2::IonImageGrabStrategyType::None;
}

m2::BaselineCorrectionType m2Data::GuiToBaselineCorrectionStrategyType()
{
  if (this->Controls()->rdBtnBaselineCorrectionTopHat->isChecked())
    return m2::BaselineCorrectionType::TopHat;
  if (this->Controls()->rdBtnBaselineCorrectionMedian->isChecked())
    return m2::BaselineCorrectionType::Median;
  return m2::BaselineCorrectionType::None;
}

void m2Data::OnResetTiling()
{
	
	auto all = AllNodes();
//	unsigned int maxWidth = 0, maxHeight = 0;
	if (all->Size() == 0)
		return;

	

	
	for (auto &e : *all)
	{
	    double initP[] = {0,0,0};
	    mitk::Point3D origin(initP);
	    mitk::Point3D prevOrigin(initP);
		if (auto *image = dynamic_cast<m2::MSImageBase *>(e->GetData()))
		{
			prevOrigin = image->GetGeometry()->GetOrigin();
			origin = image->GetGeometry()->GetOrigin();
			
			origin[0] = image->GetPropertyValue<double>("origin x");
			origin[1] = image->GetPropertyValue<double>("origin y");
			origin[2] = 0;
			image->GetGeometry()->SetOrigin(origin);
			for (auto &&kv : image->GetImageArtifacts())
				kv.second->GetGeometry()->SetOrigin(origin);
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

void m2Data::OnNextIonImage()
{
  if (m_IonImageReference)
  {
    auto mz = m_IonImageReference->mz;
    auto tol = m_Controls.spnBxTol->value();
    if (m_Controls.rbtnTolPPM->isChecked())
    {
      tol = tol * 10e-6 * mz;
    }
    this->OnGrabIonImage(mz + tol, -1);
  }
}

void m2Data::OnPrevIonImage()
{
  if (m_IonImageReference)
  {
    auto mz = m_IonImageReference->mz;
    auto tol = m_Controls.spnBxTol->value();
    if (m_Controls.rbtnTolPPM->isChecked())
    {
      tol = tol * 10e-6 * mz;
    }
    this->OnGrabIonImage(mz - tol, -1);
  }
}

void m2Data::OnProcessingNodesRequested(const QString &id)
{
  auto nodes = AllNodes();
  ApplySettingsToNodes(nodes);
  emit m2::CommunicationService::Instance()->SendProcessingNodes(id, nodes);
}

void m2Data::EmitIonImageReference()
{
  auto nodes = this->AllNodes();
  for (auto &n : *nodes)
  {
    if (auto msImage = dynamic_cast<m2::ImzMLMassSpecImage *>(n->GetData()))
    {
      if (auto current = m_IonImageReference)
      {
        const auto &mzAxis = msImage->MassAxis();
        if (current->mz >= mzAxis.front() && current->mz <= mzAxis.back())
        {
          current->name = this->m_Controls.textIonName->text().toStdString();
          msImage->GetIonImageReferenceVector().push_back(current);

          MITK_INFO << "Add IonImageReference to ion list";
          MITK_INFO << "m/z\t" << current->mz;
          MITK_INFO << "tol\t" << current->tol;
          MITK_INFO << "name\t" << current->name;
        }
      }
      else
      {
        MITK_ERROR << "Ion image ref is not set. First grab an new ion image.";
      }
    }
  }

  ApplySettingsToNodes(nodes);
  emit m2::CommunicationService::Instance()->SendProcessingNodes(QString::fromStdString(m2Data::VIEW_ID), nodes);
}

void m2Data::OnGrabIonImageFinished(mitk::DataNode * /*dataNode*/, mitk::Image * /*ionImage*/) {}

m2Data::NodesVectorType::Pointer m2Data::AllNodes()
{
  // Create a container vor all nodes that should be processed
  auto nodesToProcess = itk::VectorContainer<unsigned int, mitk::DataNode::Pointer>::New();
  if (!m_Controls.checkBoxOnlyComboBoxSelection->isChecked())
  { // is checked: all nodes that fit the predicate are processed
    auto processableNodes = GetDataStorage()->GetSubset(this->m_MassSpecDataNodeSelectionWidget->GetNodePredicate());
    if (!m_Controls.checkBoxOnlyVisible->isChecked())
    {
      for (auto node : processableNodes->CastToSTLConstContainer())
      {
        nodesToProcess->push_back(node);
      }
    }
    else
    {
      for (auto node : processableNodes->CastToSTLConstContainer())
      {
        bool visibility;
        if (node->GetBoolProperty("visible", visibility))
          if (visibility)
            nodesToProcess->push_back(node);
      }
    }
  }
  else
  { // or only the current index is used
    nodesToProcess->clear();
    auto nodes = this->m_MassSpecDataNodeSelectionWidget->GetSelectedNodes();
    if (!nodes.empty())
    {
      for (auto n : nodes)
      {
        nodesToProcess->push_back(n);
      }
    }
  }
  return nodesToProcess;
}

void m2Data::ApplySettingsToNodes(m2Data::NodesVectorType::Pointer v)
{
  for (auto dataNode : *v)
  {
    if (auto data = dynamic_cast<m2::ImzMLMassSpecImage *>(dataNode->GetData()))
      ApplySettingsToImage(data); // ImzML
                                  /*else if (auto data = dynamic_cast<mitk::Imzs*>(dataNode->GetData()))
                                    UpdateImzMLDataObject(data);*/
  }
}

void m2Data::ApplySettingsToImage(m2::ImzMLMassSpecImage *data)
{
  if (data)
  {
    using m2::ImzMLMassSpecImage;

    data->SetNormalizationStrategy(GuiToNormalizationStrategyType());
	data->SetBaselineCorrectionStrategy(GuiToBaselineCorrectionStrategyType());
    data->UseSmoothing = m_Controls.checkBoxSmoothing->isChecked();

    data->SetIonImageGrabStrategy(GuiToIonImageGrabStrategyType());
    data->SetSmoothingHalfWindowSize(m_Controls.spnBxSmoothing->value());
    data->SetBaseLinecorrectionHalfWindowSize(m_Controls.spnBxBaseline->value());

    data->SetPeakPickingBinningTolerance(m_Controls.spnBxPeakBinning->value());
  }
}

void m2Data::OnEqualizeLW()
{
  auto nodes = this->m_MassSpecDataNodeSelectionWidget->GetSelectedNodes();
  if (nodes.size() == 1)
  {
    auto node = nodes.front();
    mitk::LevelWindow lw_ref, lw_ref_normalization;
    node->GetLevelWindow(lw_ref);
    auto derivations = this->GetDataStorage()->GetDerivations(node);
    auto nPos = std::find_if(derivations->begin(), derivations->end(), [](mitk::DataNode::Pointer v) {
      return v->GetName().find("normalization") != std::string::npos;
    });
    (*nPos)->GetLevelWindow(lw_ref_normalization);

    auto allNodes = this->AllNodes();
    for (auto &&n : *allNodes)
    {
      n->SetLevelWindow(lw_ref);
      derivations = this->GetDataStorage()->GetDerivations(n);
      nPos = std::find_if(derivations->begin(), derivations->end(), [](mitk::DataNode::Pointer v) {
        return v->GetName().find("normalization") != std::string::npos;
      });

      (*nPos)->SetLevelWindow(lw_ref_normalization);
    }
  }
}

void m2Data::OnApplyTiling()
{
  auto v = m_Controls.mosaicRows->value();
  if (v < 1)
    return;
  auto all = AllNodes();
  std::vector<mitk::DataNode::Pointer> nodes;
  unsigned int maxWidth = 0, maxHeight = 0;
  if (all->Size() == 0)
    return;

  for (auto &&e : *all)
  {
    if (auto *image = dynamic_cast<m2::MSImageBase *>(e->GetData()))
    {
      maxWidth = std::max(maxWidth, image->GetDimensions()[0]);
      maxHeight = std::max(maxHeight, image->GetDimensions()[1]);
      nodes.push_back(e);
    }
  }

  int R = v;
  int N = (nodes.size() + 1) / R;
  int i = 0;
  if (N < 1)
    return;

  std::sort(nodes.begin(), nodes.end(), [](mitk::DataNode::Pointer &a, mitk::DataNode::Pointer &b) -> bool {
    return a->GetName().compare(b->GetName()) < 0;
  });

  for (auto &&e : nodes)
  {
    mitk::Point3D origin, prevOrigin;
    prevOrigin.Fill(0);
    if (auto *image = dynamic_cast<m2::MSImageBase *>(e->GetData()))
    {
      prevOrigin = image->GetGeometry()->GetOrigin();
      origin[0] = maxWidth * int(i % N) * image->GetGeometry()->GetSpacing()[0];
      origin[1] = maxHeight * int(i / N) * image->GetGeometry()->GetSpacing()[1];
      origin[2] = double(0.0);
      ++i;
      image->GetGeometry()->SetOrigin(origin);
      for (auto &&kv : image->GetImageArtifacts())
        kv.second->GetGeometry()->SetOrigin(origin);
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
  // this->RequestRenderWindowUpdate();
}

// void m2Data::OnCopyRescale()
//{
//  if (auto current = m_MassSpecDataNodeSelectionWidget->GetSelectedNode())
//  {
//    if (auto image = dynamic_cast<mitk::Image *>(current->GetData()))
//    {
//      auto node = mitk::DataNode::New();
//
//      AccessByItk(
//        image, ([&node](auto itkImage) {
//          using OutType = itk::Image<unsigned char, 3>;
//          auto rescaleFilter =
//            itk::RescaleIntensityImageFilter<typename std::remove_pointer<decltype(itkImage)>::type, OutType>::New();
//          rescaleFilter->SetInput(itkImage);
//          rescaleFilter->Update();
//          mitk::Image::Pointer out;
//          mitk::CastToMitkImage(rescaleFilter->GetOutput(), out);
//          node->SetData(out);
//        }));
//      node->SetName(current->GetName() + "_rescaled");
//      this->GetDataStorage()->Add(node);
//    }
//  }
//}

void m2Data::UpdateColorBarAndRenderWindows()
{
  mitk::ColorBarAnnotation::Pointer cbAnnotation;
  auto lookuptabel = mitk::LookupTableProperty::New();

  if (m_Controls.scaleBarEnableIon->isChecked())
  {
    cbAnnotation = m_ColorBarAnnotations[0];

    auto renderer = GetRenderWindowPart()->GetQmitkRenderWindow("axial")->GetRenderer();
    mitk::LayoutAnnotationRenderer::AddAnnotation(cbAnnotation, renderer);

    // auto processNodes = this->AllNodes();
    // for (auto &&n : *processNodes)
    //{
    //  auto derivations = this->GetDataStorage()->GetDerivations(n);
    //  for (auto &&d : *derivations)
    //  {
    //    if (d->IsVisible(nullptr))
    //    {
    //      d->SetVisibility(false);
    //    }
    //  }
    //}
    auto nodes = this->m_MassSpecDataNodeSelectionWidget->GetSelectedNodes();
    if (nodes.size() == 1)
    {
      if (nodes.front()->GetProperty(lookuptabel, "LookupTable"))
        cbAnnotation->SetLookupTable(lookuptabel->GetValue()->GetVtkLookupTable());
    }
  }
  else
  {
    cbAnnotation = m_ColorBarAnnotations[1];

    auto renderer = GetRenderWindowPart()->GetQmitkRenderWindow("axial")->GetRenderer();
    mitk::LayoutAnnotationRenderer::AddAnnotation(cbAnnotation, renderer);

    mitk::DataNode::Pointer node;
    auto processNodes = this->AllNodes();
    for (auto &&n : *processNodes)
    {
      auto derivations = this->GetDataStorage()->GetDerivations(n);
      for (auto &&n : *derivations)
      {
        if (n->GetName().find("normalization") != std::string::npos)
        {
          n->SetVisibility(true);
          node = n;
        }
      }
    }
    if (node)
    {
      if (node->GetProperty(lookuptabel, "LookupTable"))
      {
        cbAnnotation->SetLookupTable(lookuptabel->GetValue()->GetVtkLookupTable());
      }
    }
  }
}

void m2Data::OnGrabIonImage(qreal mz, qreal tol)
{
  const bool initializeNew = m2Data::Controls()->chkBxInitNew->isChecked();

  // get the selection
  auto nodesToProcess = this->AllNodes();
  if (nodesToProcess->empty())
    return;

  if (tol < 0)
  {
    tol = Controls()->spnBxTol->value();
    bool isPpm = Controls()->rbtnTolPPM->isChecked();
    tol = isPpm ? tol * 10e-6 * mz : tol;
  }
  emit m2::CommunicationService::Instance()->MassRangeChanged(mz, tol);

  auto flag = std::make_shared<std::atomic<bool>>(0);
  QString labelText = str(boost::format("%.2f +/- %.2f Da") % mz % tol).c_str();
  if (nodesToProcess->size() == 1)
    labelText = QString(nodesToProcess->front()->GetName().c_str()) + "\n" + labelText;

  this->m_Controls.labelIonImageInfo->setWordWrap(true);
  this->m_Controls.labelIonImageInfo->setText(labelText);
  this->m_Controls.spnBxMz->setValue(mz);

  this->UpdateTextAnnotations(labelText.toStdString());
  this->ApplySettingsToNodes(nodesToProcess);

  // all images share a new created ion image reference object
  // this ensures unified coloring across multiple images
  // the object is set as the currentIonImageReference during the
  // OnIonImageGrab method call. An extra action is required to add this
  // current object to the permament ion image reference vector.
  m_IonImageReference = m2::IonImageReference::New(mz, tol, "");

  for (mitk::DataNode::Pointer dataNode : *nodesToProcess)
  {
    if (m2::MSImageBase::Pointer data = dynamic_cast<m2::MSImageBase *>(dataNode->GetData()))
    {
      if (!data->IsInitialized())
        mitkThrow() << "Trying to grab an ion image but data access was not initialized properly!";

      if (dynamic_cast<m2::ImzMLMassSpecImage *>(data.GetPointer()))
      {
        if (mz > data->MassAxis().back() || mz < data->MassAxis().front())
          continue;
      }

      mitk::Image::Pointer maskImage;
      auto maskNode = FindChildNodeWithNameContaining(dataNode, "mask");

      if (maskNode && Controls()->chkBxUseMask->isChecked())
        maskImage = dynamic_cast<mitk::Image *>(maskNode->GetData());

      // if it is null, it's the first time an ion image is grabbed - disable visibility
      // so that the level window works fine for the ion image and not for the mask
      if (maskNode && maskNode->GetProperty("m/z") == nullptr)
        maskNode->SetVisibility(false);

      // The smartpointer will stay alive until all captured copies are relesed. Additional
      // all connected signals must be disconnected to make sure that the future is not kept
      // alive after the 'finished-callback' is processed.
      auto future = std::make_shared<QFutureWatcher<mitk::Image::Pointer>>();
      data->SetCurrentIonImageReference(m_IonImageReference);

      //*************** Worker Finished Callback ******************//
      // capture holds a copy of the smartpointer, so it will stay alive. Make the lambda mutable to
      // allow the manipulation of captured varaibles that are copied by '='.
      const auto futureFinished = [future, dataNode, nodesToProcess, labelText, this]() mutable {
        auto ionImage = future->result();

        if (this->m_Controls.chkBxInitNew->isChecked())
        {
          using SourceImageType = itk::Image<m2::IonImagePixelType, 3>;
          auto Caster = [&ionImage](auto *itkImage) {
            SourceImageType::Pointer srcIonImage;

            mitk::CastToItkImage(ionImage, srcIonImage);
            using ImageType = typename std::remove_pointer<decltype(itkImage)>::type;
            using FilterType = itk::RescaleIntensityImageFilter<SourceImageType, ImageType>;
            auto filter = FilterType::New();
            filter->SetInput(srcIonImage);
            filter->SetOutputMinimum(0);
            filter->SetOutputMaximum(std::numeric_limits<typename ImageType::PixelType>::max());
            filter->Update();
            mitk::CastToMitkImage(filter->GetOutput(), ionImage);
          };

          auto id = m_Controls.cmbBxType->currentIndex();
          if (id == 0)
          {
            // do nothing
          }
          else if (id == 1)
          {
            using TargetImageType = itk::Image<unsigned char, 3>;
            TargetImageType::Pointer itkIonImage;
            Caster(itkIonImage.GetPointer());
          }
          else if (id == 2)
          {
            using TargetImageType = itk::Image<unsigned short, 3>;
            TargetImageType::Pointer itkIonImage;
            Caster(itkIonImage.GetPointer());
          }
          else if (id == 3)
          {
            using TargetImageType = itk::Image<unsigned int, 3>;
            TargetImageType::Pointer itkIonImage;
            Caster(itkIonImage.GetPointer());
          }

          auto dataNode = mitk::DataNode::New();
          auto lut = mitk::LookupTable::New();
          lut->SetType(mitk::LookupTable::LookupTableType::GRAYSCALE_TRANSPARENT);
          dataNode->SetProperty("LookupTable", mitk::LookupTableProperty::New(lut));
          dataNode->SetData(ionImage);
          this->GetDataStorage()->Add(dataNode);
          if (nodesToProcess->size() == 1)
            dataNode->SetName(labelText.toStdString());
          else
            dataNode->SetName(dataNode->GetName() + "\n" + labelText.toStdString());
        }
        else
        {
          UpdateLevelWindow(dataNode);
        }
        dataNode->SetProperty("mz", ionImage->GetProperty("mz"));
        dataNode->SetProperty("tol", ionImage->GetProperty("tol"));
        mitk::LevelWindow lw;
        dataNode->GetLevelWindow(lw);
        lw.SetToImageRange(ionImage);
        dataNode->SetLevelWindow(lw);

        this->UpdateColorBarAndRenderWindows();
        this->RequestRenderWindowUpdate();

        future->disconnect();
      };

      //*************** Worker Block******************//
      const auto futureWorker = [mz, tol, data, maskImage, initializeNew, this]() {
        mitk::Timer t("Grab IonImage@[" + std::to_string(mz) + " " + std::to_string(tol) + "]");
        if (initializeNew)
        {
          auto geom = data->GetGeometry()->Clone();
          auto image = mitk::Image::New();
          image->Initialize(mitk::MakeScalarPixelType<m2::IonImagePixelType>(), *geom);
          data->GrabIonImage(mz, tol, maskImage, image);
          return image;
        }
        else
        {
          data->GrabIonImage(mz, tol, maskImage, data);
          mitk::Image::Pointer imagePtr = data.GetPointer();
          return imagePtr;
        }
      };

      //*************** Start Worker ******************//

      connect(future.get(), &QFutureWatcher<mitk::Image::Pointer>::finished, future.get(), futureFinished);
      future->setFuture(QtConcurrent::run(&m_pool, futureWorker));

      // if the future main process finished, we need to update the data nodes and their properties to
      // correctly visualize the new ion image. It requires an update on the level window propery and
      // a 'new init' if the data node is binary.

      // future->waitForFinished();
    }
  }
}

void m2Data::UpdateTextAnnotations(std::string text)
{
  static const std::array<std::string, 3> windownames = {"axial", "sagittal", "coronal"};
  if (m_TextAnnotations.size() != 3)
  {
    m_TextAnnotations.clear();
    for (int i = 0; i < 3; i++)
    {
      m_TextAnnotations.push_back(mitk::TextAnnotation2D::New());
      auto renderer = GetRenderWindowPart()->GetQmitkRenderWindow(windownames[i].c_str())->GetRenderer();
      mitk::LayoutAnnotationRenderer::AddAnnotation(m_TextAnnotations.back(), renderer);
      m_TextAnnotations.back()->SetFontSize(15);
      float color[] = {0.7, 0.7, 0.7};
      m_TextAnnotations.back()->SetFontSize(20);

      m_TextAnnotations.back()->SetColor(color);
    }
  }
  for (auto anno : m_TextAnnotations)
  {
    anno->SetText(text);
  }
}

mitk::DataNode::Pointer m2Data::FindChildNodeWithNameContaining(mitk::DataNode::Pointer &parent,
                                                                const std::string &subStr)
{
  auto deriv = this->GetDataStorage()->GetDerivations(parent.GetPointer(), nullptr, true);
  for (auto p : *deriv)
    if (p->GetName().find(subStr) != std::string::npos)
      return p;
  return nullptr;
}

void m2Data::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer> &nodes)
{
  if (!nodes.empty())
  {
    auto node = nodes.front();
    if (auto image = dynamic_cast<m2::MSImageBase *>(node->GetData()))
    {
      if (auto current = image->GetCurrentIonImageReference())
      {
        QString labelText = str(boost::format("%.2f +/- %.2f Da") % current->mz % current->tol).c_str();

        labelText += "\n";
        if (nodes.size() == 1)
          labelText += node->GetName().c_str();

        this->m_Controls.labelIonImageInfo->setWordWrap(true);
        this->m_Controls.labelIonImageInfo->setText(labelText);
        this->UpdateTextAnnotations(labelText.toStdString());
      }
    }
  }
}

void m2Data::UpdateLevelWindow(const mitk::DataNode *node)
{
  if (auto msImageBase = dynamic_cast<m2::MSImageBase *>(node->GetData()))
  {
    mitk::LevelWindow lw;
    node->GetLevelWindow(lw);
    lw.SetAuto(msImageBase);
    const_cast<mitk::DataNode *>(node)->SetLevelWindow(lw);
  }
}

void m2Data::NodeAdded(const mitk::DataNode *node)
{
  if (auto msImageBase = dynamic_cast<m2::MSImageBase *>(node->GetData()))
  {
    auto lut = mitk::LookupTable::New();
    lut->SetType(mitk::LookupTable::LookupTableType::VIRIDIS_TRANSPARENT);
    const_cast<mitk::DataNode *>(node)->SetProperty("LookupTable", mitk::LookupTableProperty::New(lut));
    const_cast<mitk::DataNode *>(node)->SetBoolProperty("binary", false);
    mitk::LevelWindow lw;
    const_cast<mitk::DataNode *>(node)->GetLevelWindow(lw);
    lw.SetToMaxWindowSize();
    const_cast<mitk::DataNode *>(node)->SetLevelWindow(lw);

    // If an MS base image object is add to the datastorage
    if (auto msImageImzML = dynamic_cast<m2::ImzMLMassSpecImage *>(msImageBase))
    {
      if (!msImageImzML->GetImageAccessInitialized())
      {
        // load associated data
        {
          auto source = msImageImzML->GetSourceList().front();
          
          auto path = source._ImzMLDataPath;

          itksys::SystemTools::ReplaceString(path, ".imzML", ".nrrd");
          if (itksys::SystemTools::FileExists(path))
          {
            source._MaskDataPath = path;
            auto data = mitk::IOUtil::Load(source._MaskDataPath).at(0);
            msImageImzML->GetImageArtifacts()["mask"] = dynamic_cast<mitk::Image *>(data.GetPointer());
            msImageImzML->PreventMaskImageInitializationOn();
          }

          path = source._ImzMLDataPath;
          itksys::SystemTools::ReplaceString(path, ".imzML", ".mps");
          if (itksys::SystemTools::FileExists(path))
          {
            source._PointsDataPath = path;
            auto data = mitk::IOUtil::Load(source._PointsDataPath).at(0);
            msImageImzML->GetImageArtifacts()["references"] = dynamic_cast<mitk::PointSet *>(data.GetPointer());
          }
        }

        // preprocessing options defined in gui passed to ImzML object
        this->ApplySettingsToImage(msImageImzML);
        // now initialize the MS image
        msImageImzML->InitializeImageAccess();

        // and add as child node
        if (msImageImzML->GetImageArtifacts().find("references") != std::end(msImageImzML->GetImageArtifacts()))
        {
          auto helperNode = mitk::DataNode::New();
          helperNode->SetName(node->GetName());
          helperNode->SetVisibility(false);
          helperNode->SetData(msImageImzML->GetImageArtifacts()["references"]);
          helperNode->SetFloatProperty("point 2D size",
                                       msImageImzML->GetGeometry()->GetSpacing()[0] * 3);              // x spacing
          helperNode->SetFloatProperty("pointsize", msImageImzML->GetGeometry()->GetSpacing()[0] * 3); // x spacing
          this->GetDataStorage()->Add(helperNode, const_cast<mitk::DataNode *>(node));
        }

        // -------------- add Mask to datastorage --------------
        {
          auto labelSetImag = mitk::ConvertImageToLabelSetImage(msImageImzML->GetMaskImage());
          msImageImzML->GetImageArtifacts()["mask"] = labelSetImag;

          auto helperNode = mitk::DataNode::New();
          helperNode->SetName(node->GetName() + "_mask");
          helperNode->SetVisibility(true);
          helperNode->SetData(msImageImzML->GetMaskImage());

          this->GetDataStorage()->Add(helperNode, const_cast<mitk::DataNode *>(node));

          labelSetImag->GetLabel(1)->SetOpacity(0);
          labelSetImag->GetActiveLabelSet()->UpdateLookupTable(1);
        }

        // -------------- add Index to datastorage --------------
        {
          auto helperNode = mitk::DataNode::New();
          helperNode->SetName(node->GetName() + "_index");
          helperNode->SetVisibility(false);
          helperNode->SetData(msImageImzML->GetIndexImage());

          this->GetDataStorage()->Add(helperNode, const_cast<mitk::DataNode *>(node));
        }

        // -------------- add Normalization to datastorage --------------
        {
          auto helperNode = mitk::DataNode::New();
          helperNode->SetName(node->GetName() + "_normalization");
          helperNode->SetVisibility(false);
          helperNode->SetData(msImageImzML->GetNormalizationImage());

          this->GetDataStorage()->Add(helperNode, const_cast<mitk::DataNode *>(node));
        }

        // 20201023: custom selection service did not work as expected
        // using NodeSelection = m2::SelectionProvider::NodeSelection;
        // using NodeWrapper = m2::SelectionProvider::NodeWrapper;
        // NodeSelection::Pointer infoSelection = NodeSelection::Pointer(new NodeSelection(NodeWrapper{node, 0}));
        // this->m_SelectionProvider->SetSelection(infoSelection);

        // 20201023: use communciation service
        emit m2::CommunicationService::Instance()->MSImageNodeAdded(node);

        // rotate 3D view
        this->RequestRenderWindowUpdate();

        if (vtkRenderWindow *renderWindow = mitk::BaseRenderer::GetRenderWindowByName("stdmulti.widget3"))
        {
          if (auto controller = mitk::BaseRenderer::GetInstance(renderWindow)->GetCameraController())
          {
            controller->SetViewToCaudal();
            /* this->GetRenderWindowPart()->GetActiveQmitkRenderWindow()->SetLayoutIndex(
               QmitkRenderWindowMenu::LayoutIndex::AXIAL);
             this->GetRenderWindowPart()->GetActiveQmitkRenderWindow()->LayoutDesignChanged(
               QmitkRenderWindowMenu::LayoutDesign::ONE_BIG);*/
          }
        }
        // this->GetRenderWindowPart()->GetActiveQmitkRenderWindow()->LayoutDesignChanged(QmitkRenderWindowMenu::LayoutDesign::ONE_LEFT_3D_RIGHT);

        // controller->SetViewToCaudal();
        // controller->SetViewToCranial();
        // controller->SetViewToAnterior();

        //      controller->SetViewToPosterior();
        //      controller->SetViewToSinister();
      }
      // -------------- end read operation --------------
    }
  }
}

void m2Data::NodeRemoved(const mitk::DataNode *node)
{
  if (auto msImageBase = dynamic_cast<m2::MSImageBase *>(node->GetData()))
  {
    if (dynamic_cast<m2::ImzMLMassSpecImage *>(msImageBase))
    {
      auto derivations = this->GetDataStorage()->GetDerivations(node);
      for (auto &&d : *derivations)
      {
        this->GetDataStorage()->Remove(d);
      }
    }
  }
}
