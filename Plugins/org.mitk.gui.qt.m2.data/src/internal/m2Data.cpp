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

#include "Qm2OpenSlideImageIOHelperDialog.h"
#include <QColorDialog>
#include <QComboBox>
#include <QInputDialog>
#include <QmitkRenderWindow.h>
#include <QtConcurrent>
#include <boost/format.hpp>
#include <itkRescaleIntensityImageFilter.h>
#include <itksys/SystemTools.hxx>
#include <m2FsmSpectrumImage.h>
#include <m2ImzMLSpectrumImage.h>
#include <m2IntervalVector.h>
#include <m2SpectrumImage.h>
#include <m2SpectrumImageDataInteractor.h>
#include <m2Process.hpp>
#include <m2SpectrumImageStack.h>
#include <m2ShiftMapImageFilter.h>
#include <m2SubdivideImage2DFilter.h>
#include <m2UIUtils.h>
#include <mitkColorProperty.h>
#include <mitkCoreServices.h>
#include <mitkIPreferences.h>
#include <mitkIPreferencesService.h>
#include <mitkImageCast.h>
#include <mitkIOUtil.h>
#include <mitkLayoutAnnotationRenderer.h>
#include <mitkLookupTableProperty.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateOr.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateProperty.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkImageVtkMapper2D.h>
#include <mitkImageAccessByItk.h>
#include <regex>
// #include <Qm2AssociatedFilesDialog.h>

const std::string m2Data::VIEW_ID = "org.mitk.views.m2.data";

void m2Data::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi(parent);
  m_Parent = parent;

  InitBaselineCorrectionControls();
  InitNormalizationControls();
  InitRangePoolingControls();
  InitSmoothingControls();
  InitIntensityTransformationControls();
  InitImageNormalizationControls();
  InitImageSmoothingControls();

  auto serviceRef = m2::UIUtils::Instance();
  connect(serviceRef, SIGNAL(UpdateImage(qreal, qreal)), this, SLOT(OnGenerateImageData(qreal, qreal)));

  connect(m_Controls.btnCreateImage,
          &QAbstractButton::clicked,
          this,
          [&] { OnGenerateImageData(m_Controls.spnBxMz->value(), FROM_GUI); });

  // step through
  // signals are triggered by key strokes (arrows) from spectrum view (Spectrum.ccp)
  auto UIUtilsObject = m2::UIUtils::Instance();
  connect(UIUtilsObject, SIGNAL(PreviousImage()), this, SLOT(OnCreatePrevImage()));
  connect(UIUtilsObject, SIGNAL(NextImage()), this, SLOT(OnCreateNextImage()));
  connect(UIUtilsObject, SIGNAL(PreviousPeakImage()), this, SLOT(OnCreatePrevPeakImage()));
  connect(UIUtilsObject, SIGNAL(NextPeakImage()), this, SLOT(OnCreateNextPeakImage()));
  connect(UIUtilsObject, SIGNAL(IncreaseTolerance()), this, SLOT(OnIncreaseTolerance()));
  connect(UIUtilsObject, SIGNAL(DecreaseTolerance()), this, SLOT(OnDecreaseTolerance()));

  connect(m_Controls.btnCreateShiftMap, SIGNAL(clicked()), this, SLOT(OnCreateShiftMap()));


  // Position list/history
  // connect(m_Controls.listWidgetPositions, &QTableWidget::cellDoubleClicked, this, [this](int row, int){
  //   auto text = m_Controls.listWidgetPositions->item(row, 0)->text();
  //   if (text == "Empty")
  //     return;
  //   if (row == 0){
  //     m_Controls.listWidgetPositions->insertRow(0);
  //     m_Controls.listWidgetPositions->setItem(0, 0, new QTableWidgetItem("Empty")); 
  //   }
  //   else
  //   {
  //     auto x = text.split(" ")[1].toDouble();
  //     auto tol = text.split(" ")[3].toDouble();
  //     this->OnGenerateImageData(x, tol);
  //   }
  // });
  

  // Settings (show hlper objects)
  using namespace mitk;

  const auto toggleByType = [&](bool isChecked, m2::NormalizationStrategyType type)
  {
    const auto a = TNodePredicateDataType<m2::SpectrumImageStack>::New();
    const auto b = TNodePredicateDataType<m2::ImzMLSpectrumImage>::New();
    const auto c = TNodePredicateDataType<m2::FsmSpectrumImage>::New();
    const auto predicate = NodePredicateOr::New(a, b, c);
    const auto nodes = this->GetDataStorage()->GetSubset(predicate);
    for (const auto &node : *nodes)
    {
      if (auto image = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
      {
        if (!image->GetNormalizationImageStatus(type))
          image->InitializeNormalizationImage(type);

        std::string name = "NormalizationImage" + m2::to_string(type);

        const auto derivations = this->GetDataStorage()->GetDerivations(node);
        for (const auto &dNode : *derivations)
        {
          if (dNode->GetName().find(name) != std::string::npos)
          {
            dNode->SetVisibility(isChecked);

            if (!isChecked)// to keep the data manger clean hide the DataNode by tagging it as helper object
              dNode->SetBoolProperty("helper object", true);
            else
            {
              // if it should be visible in the data storage remove the helper object property
              dNode->RemoveProperty("helper object");

              // update level window
              mitk::LevelWindow lw;
              dNode->GetLevelWindow(lw);
              lw.SetAuto(image->GetNormalizationImage(type));
              const auto max = lw.GetRangeMax();
              lw.SetWindowBounds(0, max);
              dNode->SetLevelWindow(lw);
            }
            this->RequestRenderWindowUpdate();
          }
        }
      }
    }
  };

  const auto toggleByName = [&](bool isChecked, const char *name)
  {
    const auto a = TNodePredicateDataType<m2::SpectrumImageStack>::New();
    const auto b = TNodePredicateDataType<m2::ImzMLSpectrumImage>::New();
    const auto c = TNodePredicateDataType<m2::FsmSpectrumImage>::New();
    const auto predicate = NodePredicateOr::New(a, b, c);
    const auto nodes = this->GetDataStorage()->GetSubset(predicate);
    for (const auto &node : *nodes)
    {
      const auto derivations = this->GetDataStorage()->GetDerivations(node);
      for (const auto &dNode : *derivations)
      {
        if (dNode->GetName().find(name) != std::string::npos)
        {
          dNode->SetVisibility(isChecked);
          if (!isChecked)
            dNode->SetBoolProperty("helper object", true);
          else
            dNode->RemoveProperty("helper object");
          this->RequestRenderWindowUpdate();
        }
      }
    }
  };

  // connect checkboxes to handle visibility of helper objects in the DataManager
  int i = 1;
  for (m2::NormalizationStrategyType type : m2::NormalizationStrategyTypeList)
  {
    auto ckBox = new QCheckBox(("Show " + m2::to_string(type) + " normalization images").c_str(), m_Controls.settings);
    QHBoxLayout *layout = (QHBoxLayout *)(m_Controls.settings->layout());
    layout->insertWidget(layout->indexOf(m_Controls.hLineNormImages) + i, ckBox);
    ckBox->setObjectName(("ckBoxNormalizationImage" + m2::to_string(type)).c_str());

    connect(ckBox, &QCheckBox::toggled, this, [toggleByType, type](int v) { toggleByType(v, type); });
    ++i;
  }


  connect(m_Controls.showIndexImages,
          &QCheckBox::toggled,
          this,
          [toggleByName](bool isChecked) { toggleByName(isChecked, "IndexImage"); });
  connect(m_Controls.showMaskImages,
          &QCheckBox::toggled,
          this,
          [toggleByName](bool isChecked) { toggleByName(isChecked, "MaskImage"); });
  connect(m_Controls.showMeanSpectrum,
          &QCheckBox::toggled,
          this,
          [toggleByName](bool isChecked) { toggleByName(isChecked, "MeanSpectrum"); });
  connect(m_Controls.showMaxSpectrum,
          &QCheckBox::toggled,
          this,
          [toggleByName](bool isChecked) { toggleByName(isChecked, "MaxSpectrum"); });
  connect(m_Controls.showSingleSpectrum,
          &QCheckBox::stateChanged,
          this,
          [toggleByName](bool isChecked) { toggleByName(isChecked, "SingleSpectrum"); });
  connect(m_Controls.showCentroidSpectrum,
          &QCheckBox::stateChanged,
          this,
          [toggleByName](bool isChecked) { toggleByName(isChecked, "CentroidSpectrum"); });

  connect(m2::UIUtils::Instance(),
          &m2::UIUtils::RequestTolerance,
          this,
          [this](float x, float &tol)
          {
            tol = m_Controls.spnBxTol->value();
            if (m_Controls.rbtnTolPPM->isChecked())
              tol = m2::PartPerMillionToFactor(tol) * x;
          });

  // Imaging controls
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();

  connect(m_Controls.spnBxTol,
          qOverload<double>(&QDoubleSpinBox::valueChanged),
          this,
          [this, preferences](int)
          {
            auto value = m_Controls.spnBxTol->value();
            preferences->PutFloat("m2aia.signal.Tolerance", value);
          });

  connect(m_Controls.CBNormalization,
          qOverload<int>(&QComboBox::currentIndexChanged),
          this,
          [this, preferences](int)
          {
            auto value = m_Controls.CBNormalization->currentData().toUInt();
            preferences->PutInt("m2aia.signal.NormalizationStrategy", value);
          });

  connect(m_Controls.CBTransformation,
          qOverload<int>(&QComboBox::currentIndexChanged),
          this,
          [this, preferences](int)
          {
            auto value = m_Controls.CBTransformation->currentData().toUInt();
            preferences->PutInt("m2aia.signal.IntensityTransformationStrategy", value);
          });

  connect(m_Controls.CBImagingStrategy,
          qOverload<int>(&QComboBox::currentIndexChanged),
          this,
          [this, preferences](int)
          {
            auto value = m_Controls.CBImagingStrategy->currentData().toUInt();
            preferences->PutInt("m2aia.signal.RangePoolingStrategy", value);
          });
  connect(m_Controls.CBSmoothing,
          qOverload<int>(&QComboBox::currentIndexChanged),
          this,
          [this, preferences](int)
          {
            auto value = m_Controls.CBSmoothing->currentData().toUInt();
            preferences->PutInt("m2aia.signal.SmoothingStrategy", value);
          });
  connect(m_Controls.CBImageNormalization,
          qOverload<int>(&QComboBox::currentIndexChanged),
          this,
          [this, preferences](int)
          {
            auto value = m_Controls.CBImageNormalization->currentData().toUInt();
            preferences->PutInt("m2aia.signal.ImageNormalizationStrategy", value);
          });
  connect(m_Controls.CBImageSmoothing,
          qOverload<int>(&QComboBox::currentIndexChanged),
          this,
          [this, preferences](int)
          {
            auto value = m_Controls.CBImageSmoothing->currentData().toUInt();
            preferences->PutInt("m2aia.signal.ImageSmoothingStrategy", value);
          });


          

  // default values
  m_Controls.spnBxTol->setValue(preferences->GetFloat("m2aia.signal.Tolerance", 75));
  
  m_Controls.CBNormalization->setCurrentIndex(
    preferences->GetInt("m2aia.signal.NormalizationStrategy", to_underlying(m2::NormalizationStrategyType::None)));
  
  m_Controls.CBTransformation->setCurrentIndex(preferences->GetInt("m2aia.signal.IntensityTransformationStrategy",
                                                                   to_underlying(m2::NormalizationStrategyType::None)));
  m_Controls.CBSmoothing->setCurrentIndex(
    preferences->GetInt("m2aia.signal.SmoothingStrategy", to_underlying(m2::NormalizationStrategyType::None)));
  
  m_Controls.CBImagingStrategy->setCurrentIndex(
    preferences->GetInt("m2aia.signal.RangePoolingStrategy", to_underlying(m2::NormalizationStrategyType::Mean)));
  
  m_Controls.CBImageNormalization->setCurrentIndex(
    preferences->GetInt("m2aia.signal.ImageNormalizationStrategy", to_underlying(m2::ImageNormalizationStrategyType::None)));
  
  m_Controls.CBImageNormalization->setCurrentIndex(
    preferences->GetInt("m2aia.signal.ImageSmoothingStrategy", to_underlying(m2::ImageSmoothingStrategyType::None)));

  // Make sure, that data nodes added before this view
  // is initialized are handled correctly!!
  auto nodes = this->GetDataStorage()->GetAll();
  unsigned int count = 0;
  for (auto n : *nodes)
  {
    if (auto I = dynamic_cast<m2::ImzMLSpectrumImage *>(n->GetData()))
    {
      if (I->GetImageAccessInitialized())
      {
        n->SetVisibility(false);
      }
      else
      {
        count++;
        n->SetStringProperty("UpdateRequired", "uninitialized imzML");
      }
    }
  }

  if (count)
  {
    auto future = QtConcurrent::run([parent, this]() { QThread::currentThread()->sleep(1); });

    auto watcher = std::make_shared<QFutureWatcher<void>>();
    watcher->setFuture(future);
    connect(watcher.get(),
            &QFutureWatcher<void>::finished,
            [watcher, this]() mutable
            {
              auto nodes = this->GetDataStorage()->GetAll();

              for (auto n : *nodes)
              {
                if (n->GetProperty("UpdateRequired"))
                {
                  auto I = dynamic_cast<m2::SpectrumImage *>(n->GetData());
                  this->ApplySettingsToImage(I);
                  this->GetDataStorage()->Remove(n);
                  this->GetDataStorage()->Add(n);
                }
              }
              mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
              this->RequestRenderWindowUpdate();
              watcher->disconnect();
              watcher.reset();
            });
  }


  connect(&m_ResetPreventDataStorageOverload, &QFutureWatcher<void>::finished, this, [this]() {
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
    this->RequestRenderWindowUpdate();
  });
}

void m2Data::InitToleranceControls()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto defaultValue = preferences->GetFloat("m2aia.signal.Tolerance", 75.0);
  m_Controls.spnBxTol->setValue(defaultValue);
}

void m2Data::InitNormalizationControls()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto defaultValue =
    preferences->GetInt("m2aia.signal.NormalizationStrategy", to_underlying(m2::NormalizationStrategyType::None));
  auto cb = Controls()->CBNormalization;
  for (unsigned int i = 0; i < m2::NormalizationStrategyTypeNames.size(); ++i)
    cb->addItem(m2::NormalizationStrategyTypeNames[i].c_str(), {i});
  cb->setCurrentIndex(defaultValue);
}

m2::NormalizationStrategyType m2Data::GuiToNormalizationStrategyType()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto value =
    preferences->GetInt("m2aia.signal.NormalizationStrategy", to_underlying(m2::NormalizationStrategyType::None));
  return static_cast<m2::NormalizationStrategyType>(value);
}

void m2Data::InitImageNormalizationControls(){
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto defaultValue =
    preferences->GetInt("m2aia.signal.ImageNormalizationStrategy", to_underlying(m2::ImageNormalizationStrategyType::None));
  auto cb = Controls()->CBImageNormalization;
  for (unsigned int i = 0; i < m2::ImageNormalizationStrategyTypeNames.size(); ++i)
    cb->addItem(m2::ImageNormalizationStrategyTypeNames[i].c_str(), {i});
  cb->setCurrentIndex(defaultValue);
}


m2::ImageNormalizationStrategyType m2Data::GuiToImageNormalizationStrategyType()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto value =
    preferences->GetInt("m2aia.signal.ImageNormalizationStrategy", to_underlying(m2::ImageNormalizationStrategyType::None));
  return static_cast<m2::ImageNormalizationStrategyType>(value);
}

m2::ImageSmoothingStrategyType m2Data::GuiToImageSmoothingStrategyType()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto value =
    preferences->GetInt("m2aia.signal.ImageSmoothingStrategy", to_underlying(m2::ImageSmoothingStrategyType::None));
  return static_cast<m2::ImageSmoothingStrategyType>(value);
}

void m2Data::InitImageSmoothingControls()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto defaultValue =
    preferences->GetInt("m2aia.signal.ImageSmoothingStrategy", to_underlying(m2::ImageSmoothingStrategyType::None));
  auto cb = Controls()->CBImageSmoothing;
  for (unsigned int i = 0; i < m2::ImageSmoothingStrategyTypeNames.size(); ++i)
    cb->addItem(m2::ImageSmoothingStrategyTypeNames[i].c_str(), {i});
  cb->setCurrentIndex(defaultValue);
}


void m2Data::InitIntensityTransformationControls()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto defaultValue = preferences->GetInt("m2aia.signal.IntensityTransformationStrategy",
                                          to_underlying(m2::IntensityTransformationType::None));
  // Combo Box
  auto cb = Controls()->CBTransformation;
  for (unsigned int i = 0; i < m2::IntensityTransformationTypeNames.size(); ++i)
    cb->addItem(m2::IntensityTransformationTypeNames[i].c_str(), {i});

  cb->setCurrentIndex(defaultValue);
}

m2::IntensityTransformationType m2Data::GuiToIntensityTransformationStrategyType()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto value = preferences->GetInt("m2aia.signal.IntensityTransformationStrategy",
                                   to_underlying(m2::IntensityTransformationType::None));
  return static_cast<m2::IntensityTransformationType>(value);
}

void m2Data::InitRangePoolingControls()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto defaultValue =
    preferences->GetInt("m2aia.signal.RangePoolingStrategy", to_underlying(m2::RangePoolingStrategyType::Mean));
  auto cb = Controls()->CBImagingStrategy;
  for (unsigned int i = 0; i < m2::RangePoolingStrategyTypeNames.size(); ++i)
    cb->addItem(m2::RangePoolingStrategyTypeNames[i].c_str(), {i}); // add i as data

  cb->setCurrentIndex(defaultValue);
}

m2::RangePoolingStrategyType m2Data::GuiToRangePoolingStrategyType()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto value =
    preferences->GetInt("m2aia.signal.RangePoolingStrategy", to_underlying(m2::RangePoolingStrategyType::None));
  return static_cast<m2::RangePoolingStrategyType>(value);
}

void m2Data::InitSmoothingControls()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto defaultValue = preferences->GetInt("m2aia.signal.SmoothingStrategy", to_underlying(m2::SmoothingType::None));
  auto cb = Controls()->CBSmoothing;
  for (unsigned int i = 0; i < m2::SmoothingTypeNames.size(); ++i)
    cb->addItem(m2::SmoothingTypeNames[i].c_str(), {i});
  cb->setCurrentIndex(defaultValue);

  // Spin Box
  m_Controls.spnBxSmoothing->setValue(preferences->GetInt("m2aia.signal.SmoothingValue", 2));
  connect(m_Controls.spnBxSmoothing,
          qOverload<int>(&QSpinBox::valueChanged),
          this,
          [this, preferences](int)
          { preferences->PutInt("m2aia.signal.SmoothingValue", m_Controls.spnBxSmoothing->value()); });
}

m2::SmoothingType m2Data::GuiToSmoothingStrategyType()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto value = preferences->GetInt("m2aia.signal.SmoothingStrategy", to_underlying(m2::SmoothingType::None));
  return static_cast<m2::SmoothingType>(value);
}

void m2Data::InitBaselineCorrectionControls()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto defaultValue =
    preferences->GetInt("m2aia.signal.BaselineCorrectionStrategy", to_underlying(m2::BaselineCorrectionType::None));
  auto cb = Controls()->CBBaselineCorrection;
  for (unsigned int i = 0; i < m2::BaselineCorrectionTypeNames.size(); ++i)
    cb->addItem(m2::BaselineCorrectionTypeNames[i].c_str(), {i});
  cb->setCurrentIndex(defaultValue);

  connect(m_Controls.CBBaselineCorrection,
          qOverload<int>(&QComboBox::currentIndexChanged),
          this,
          [this, preferences](int)
          {
            auto value = m_Controls.CBBaselineCorrection->currentData().toUInt();
            preferences->PutInt("m2aia.signal.BaselineCorrectionStrategy", value);
          });

  // Spin Box
  m_Controls.spnBxBaseline->setValue(preferences->GetInt("m2aia.signal.BaselineCorrectionValue", 50));
  connect(m_Controls.spnBxBaseline,
          qOverload<int>(&QSpinBox::valueChanged),
          this,
          [this, preferences](int)
          { preferences->PutInt("m2aia.signal.BaselineCorrectionValue", m_Controls.spnBxBaseline->value()); });
}

m2::BaselineCorrectionType m2Data::GuiToBaselineCorrectionStrategyType()
{
  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
  auto *preferences = preferencesService->GetSystemPreferences();
  auto value =
    preferences->GetInt("m2aia.signal.BaselineCorrectionStrategy", to_underlying(m2::BaselineCorrectionType::None));
  return static_cast<m2::BaselineCorrectionType>(value);
}

void m2Data::OnDecreaseTolerance()
{
  m_Controls.spnBxTol->setValue(m_Controls.spnBxTol->value() * 0.9);
  OnGenerateImageData(m_Controls.spnBxMz->value(), FROM_GUI);
}

void m2Data::OnIncreaseTolerance()
{
  m_Controls.spnBxTol->setValue(m_Controls.spnBxTol->value() * 1.1);
  OnGenerateImageData(m_Controls.spnBxMz->value(), FROM_GUI);
}

void m2Data::OnCreateNextImage()
{
  auto center = m_Controls.spnBxMz->value();
  auto offset = m_Controls.spnBxTol->value();
  if (m_Controls.rbtnTolPPM->isChecked())
  {
    offset = m2::PartPerMillionToFactor(offset)*.5 * center;
  }
  this->OnGenerateImageData(center + offset, FROM_GUI);
}

void m2Data::OnCreatePrevImage()
{
  auto center = m_Controls.spnBxMz->value();
  auto offset = m_Controls.spnBxTol->value();
  if (m_Controls.rbtnTolPPM->isChecked())
  {
    offset = m2::PartPerMillionToFactor(offset)*.5 * center;
  }
  this->OnGenerateImageData(center - offset, FROM_GUI);
}

void m2Data::OnCreateNextPeakImage()
{
  auto predicate = mitk::TNodePredicateDataType<m2::IntervalVector>::New();
  auto processableNodes = GetDataStorage()->GetSubset(predicate)->CastToSTLConstContainer();

  auto center = m_Controls.spnBxMz->value();
  auto tolerance = m_Controls.spnBxTol->value();
  if (m_Controls.rbtnTolPPM->isChecked())
    tolerance = m2::PartPerMillionToFactor(tolerance)*.5 * center;
  
  std::vector<m2::Interval> nearestValues;
  for (auto node : processableNodes)
  {
    if (auto intervalVector = dynamic_cast<m2::IntervalVector *>(node->GetData()))
    {
      const auto intervalType = to_underlying(intervalVector->GetType());
      const auto centroidType = to_underlying(m2::SpectrumFormat::Centroid);
      if (intervalType & centroidType)
      {
        auto intervals = intervalVector->GetIntervals();
        auto nearestElement = std::min_element(intervals.begin(),
                                               intervals.end(),
                                               [center, tolerance](const m2::Interval &a, const m2::Interval &b)
                                               {
                                                 // Ensure both are greater than target, and compare only those
                                                 if (a.x.mean() <= center + tolerance)
                                                   return false;
                                                 if (b.x.mean() <= center + tolerance)
                                                   return true;
                                                 return a.x.mean() < b.x.mean();
                                               });
        if (nearestElement == intervals.end() || nearestElement->x.mean() <= center)
          continue;
        
        nearestValues.push_back(*nearestElement);
      }
    }
  }

  auto nearestElement = std::min_element(nearestValues.begin(),
                   nearestValues.end(),
                   [center](const m2::Interval &a, const m2::Interval &b)
                   {
                     return a.x.mean() < b.x.mean();
                   });

  if (nearestElement == nearestValues.end())
    return; // no peak found

  this->OnGenerateImageData(nearestElement->x.mean(), FROM_GUI);
}

void m2Data::OnCreatePrevPeakImage()
{
  auto predicate = mitk::TNodePredicateDataType<m2::IntervalVector>::New();
  auto processableNodes = GetDataStorage()->GetSubset(predicate)->CastToSTLConstContainer();

  auto center = m_Controls.spnBxMz->value();
  auto tolerance = m_Controls.spnBxTol->value();
  if (m_Controls.rbtnTolPPM->isChecked())
    tolerance = m2::PartPerMillionToFactor(tolerance)*.5 * center;
  
  std::vector<m2::Interval> nearestValues;
  for (auto node : processableNodes)
  {
    if (auto intervalVector = dynamic_cast<m2::IntervalVector *>(node->GetData()))
    {
      const auto intervalType = to_underlying(intervalVector->GetType());
      const auto centroidType = to_underlying(m2::SpectrumFormat::Centroid);
      if (intervalType & centroidType)
      {
        auto intervals = intervalVector->GetIntervals();
        auto nearestElement = std::min_element(intervals.begin(),
                                               intervals.end(),
                                               [center, tolerance](const m2::Interval &a, const m2::Interval &b)
                                               {
                                                 // Ensure both are greater than target, and compare only those
                                                 if (a.x.mean() >= center - tolerance)
                                                   return false;
                                                 if (b.x.mean() >= center - tolerance)
                                                   return true;
                                                 return a.x.mean() > b.x.mean();
                                               });
        if (nearestElement == intervals.end() || nearestElement->x.mean() >= center)
          continue;
        
        nearestValues.push_back(*nearestElement);
      }
    }
  }

  auto nearestElement = std::min_element(nearestValues.begin(),
                   nearestValues.end(),
                   [center](const m2::Interval &a, const m2::Interval &b)
                   {
                     return a.x.mean() > b.x.mean();
                   });

  if (nearestElement == nearestValues.end()) 
    return; // no peak found

  this->OnGenerateImageData(nearestElement->x.mean(), FROM_GUI);
}

void m2Data::ApplySettingsToNodes(m2::UIUtils::NodesVectorType::Pointer v)
{
  for (auto dataNode : *v)
  {
    if (auto data = dynamic_cast<m2::SpectrumImage *>(dataNode->GetData()))
      ApplySettingsToImage(data);
  }
}

void m2Data::ApplySettingsToImage(m2::SpectrumImage *data)
{
  if (data)
  {
    data->SetNormalizationStrategy(GuiToNormalizationStrategyType());
    data->SetBaselineCorrectionStrategy(GuiToBaselineCorrectionStrategyType());
    data->SetSmoothingStrategy(GuiToSmoothingStrategyType());
    data->SetRangePoolingStrategy(GuiToRangePoolingStrategyType());
    data->SetIntensityTransformationStrategy(GuiToIntensityTransformationStrategyType());

    data->SetImageNormalizationStrategy(GuiToImageNormalizationStrategyType());
    data->SetImageSmoothingStrategy(GuiToImageSmoothingStrategyType());

    data->SetSmoothingHalfWindowSize(m_Controls.spnBxSmoothing->value());
    data->SetBaseLineCorrectionHalfWindowSize(m_Controls.spnBxBaseline->value());
    data->SetUseToleranceInPPM(m_Controls.rbtnTolPPM->isChecked());

    // data->SetBinningTolerance(m_Controls.spnBxPeakBinning->value());
  }
}

void m2Data::OnGenerateImageData(mitk::DataNode::Pointer node,
                                 qreal xRangeCenter,
                                 qreal xRangeTol,
                                 bool emitRangeChanged)
{
  // tol < 0 indicates "use gui tol"
  if (xRangeTol < 0)
  {
    xRangeTol = Controls()->spnBxTol->value();
    bool isPpm = Controls()->rbtnTolPPM->isChecked();
    xRangeTol = isPpm ? m2::PartPerMillionToFactor(xRangeTol) * xRangeCenter : xRangeTol;
  }

  if (emitRangeChanged)
    emit m2::UIUtils::Instance()->RangeChanged(xRangeCenter, xRangeTol);

  this->m_Controls.spnBxMz->setValue(xRangeCenter);

  if (m2::SpectrumImage::Pointer data = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
  {
    auto xMin = data->GetPropertyValue<double>("m2aia.xs.min");
    auto xMax = data->GetPropertyValue<double>("m2aia.xs.max");
    if (xRangeCenter > xMax || xRangeCenter < xMin){
      mitk::ImagePixelWriteAccessor<m2::DisplayImagePixelType, 3> acc(data);
      auto N = std::accumulate(data->GetDimensions(), data->GetDimensions()+3, 1, std::multiplies<int>());
      std::fill(acc.GetData(), acc.GetData() + N, 0);
      UpdateLevelWindow(node);
      this->RequestRenderWindowUpdate();

      return;
    }

    ApplySettingsToImage(data);
    if (!data->IsInitialized())
      mitkThrow() << "Trying to grab an ion image but data access was not initialized properly!";

    mitk::Image::Pointer maskImage = data->GetMaskImage();

    // The smartpointer will stay alive until all captured copies are relesed. Additional
    // all connected signals must be disconnected to make sure that the future is not kept
    // alive after the 'finished-callback' is processed.
    auto future = std::make_shared<QFutureWatcher<mitk::Image::Pointer>>();

    //*************** Worker Finished Callback ******************//
    // capture holds a copy of the smartpointer, so it will stay alive. Make the lambda mutable to
    // allow the manipulation of captured varaibles that are copied by '='.
    const auto futureFinished = [future, node, this]() mutable
    {
      auto image = future->result();
      UpdateLevelWindow(node);
      // UpdateSpectrumImageTable(node);
      node->SetProperty("m2aia.xs.selection.center", image->GetProperty("m2aia.xs.selection.center"));
      node->SetProperty("m2aia.xs.selection.tolerance", image->GetProperty("m2aia.xs.selection.tolerance"));
      this->RequestRenderWindowUpdate();
      future->disconnect();
    };

    //*************** Worker Block******************//
    const auto futureWorker = [xRangeCenter, xRangeTol, data, maskImage, this]()
    {
      // m2::Timer t("Create image @[" + std::to_string(xRangeCenter) + " " + std::to_string(xRangeTol) + "]");
      if (m_InitializeNewNode)
      {
        auto geom = data->GetGeometry()->Clone();
        auto image = mitk::Image::New();
        image->Initialize(mitk::MakeScalarPixelType<m2::DisplayImagePixelType>(), *geom);
        data->GetImage(xRangeCenter, xRangeTol, maskImage, image);
        return image;
      }
      else
      {
        data->GetImage(xRangeCenter, xRangeTol, maskImage, data);
        mitk::Image::Pointer imagePtr = data.GetPointer();
        return imagePtr;
      }
    };

    //*************** Start Worker ******************//

    connect(future.get(), &QFutureWatcher<mitk::Image::Pointer>::finished, future.get(), futureFinished);
    future->setFuture(QtConcurrent::run(&m_pool, futureWorker));
  }
}


void m2Data::OnCreateShiftMap()
{
// get the selection
  auto nodesToProcess = m2::UIUtils::AllNodes(GetDataStorage());

  if (nodesToProcess->empty())
    return;

  // process all nodes
  for (mitk::DataNode::Pointer dataNode : *nodesToProcess)
    if (m2::ImzMLSpectrumImage::Pointer data = dynamic_cast<m2::ImzMLSpectrumImage *>(dataNode->GetData()))
    {

      auto shiftMapFilter = m2::ShiftMapImageFilter::New();

      shiftMapFilter->SetInput(data);
      shiftMapFilter->GenerateData();
      

      auto node = mitk::DataNode::New();
      node->SetData(shiftMapFilter->GetOutput(0));
      node->SetName("Absolute_mz_shift");
      GetDataStorage()->Add(node);

      node = mitk::DataNode::New();
      node->SetData(shiftMapFilter->GetOutput(1));
      node->SetName("Index_mz_shift");
      GetDataStorage()->Add(node);

    }
}


void m2Data::OnGenerateImageData(qreal xRangeCenter, qreal xRangeTol)
{
  // get the selection
  auto nodesToProcess = m2::UIUtils::AllNodes(GetDataStorage());


  if (nodesToProcess->empty())
    return;

  if (xRangeTol < 0){
    xRangeTol = Controls()->spnBxTol->value();
    bool isPpm = Controls()->rbtnTolPPM->isChecked();
    xRangeTol = isPpm ? m2::PartPerMillionToFactor(xRangeTol) * xRangeCenter : xRangeTol;
  }

  // Add element to position list if not already present
  // auto text = "m/z " + QString::number(xRangeCenter) + " ± " + QString::number(xRangeTol, 103, 3) + " Da";
  // if(m_Controls.listWidgetPositions->rowCount() == 0){
  //   m_Controls.listWidgetPositions->insertRow(0);
  //   auto item = new QTableWidgetItem(text);
  //   m_Controls.listWidgetPositions->setItem(0, 0, item);
    


  // }else{
  //   m_Controls.listWidgetPositions->itemAt(0, 0)->setText(text);
  // }


  
  emit m2::UIUtils::Instance()->RangeChanged(xRangeCenter, xRangeTol);
  this->m_Controls.spnBxMz->setValue(xRangeCenter);
  auto flag = std::make_shared<std::atomic<bool>>(0);

  QString labelText = str(boost::format("%.4f +/- %.2f Da") % xRangeCenter % xRangeTol).c_str();

  if (nodesToProcess->size() == 1)
  {
    auto node = nodesToProcess->front();
    if (auto image = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
    {
      std::string xLabel = image->GetSpectrumType().XAxisLabel;
      labelText = "[" + QString(xLabel.c_str()) + "]" + labelText;
    }
    labelText = QString(node->GetName().c_str()) + "\n" + labelText;
  }
  
  this->UpdateTextAnnotations(labelText.toStdString());

  // set GUI settings to the selected nodes
  ApplySettingsToNodes(nodesToProcess);

  // process all nodes
  for (mitk::DataNode::Pointer dataNode : *nodesToProcess)
    if (m2::SpectrumImage::Pointer data = dynamic_cast<m2::SpectrumImage *>(dataNode->GetData()))
      OnGenerateImageData(dataNode, xRangeCenter, xRangeTol, false); // do not emit
}

void m2Data::UpdateSpectrumImageTable(const mitk::DataNode *node)
{
  if (dynamic_cast<m2::SpectrumImage *>(node->GetData()))
  {
    // auto widgets = m_Parent->findChildren<Qm2ImageColorWidget *>("Qm2ImageColorWidget");
    // MITK_INFO << widgets.count();
    //  if(!widgets.count() || widgets.front()->visibilityCheckBox()->isChecked()){
    //    MITK_INFO << widgets.back();
    //    auto c = m_Controls.verticalLayout->count();
    //    m_Controls.verticalLayout->insertWidget(c-1-widgets.count(), new Qm2ImageColorWidget(m_Parent));

    //  }
  }
  //   auto item = m_Controls.tableWidget->item(0, 1);
  //   std::string xLabel = image->GetPropertyValue<std::string>("x_label");
  //   double center = image->GetPropertyValue<double>("m2aia.xs.selection.center");
  //   double tol = image->GetPropertyValue<double>("m2aia.xs.selection.tolerance");
  //   QString labelText = str(boost::format(xLabel + " %.2f +/- %.2f Da") % center % tol).c_str();
  //   item->setText(labelText);

  // }

  // item = m_Controls.tableWidget->item(0, 0)->setData()
}

void m2Data::OnRenderSpectrumImages(double min, double max)
{
  Q_UNUSED(min);
  Q_UNUSED(max);
}

void m2Data::UpdateTextAnnotations(std::string /*text*/)
{
  // static const std::array<std::string, 3> windownames = {"axial", "sagittal", "coronal"};
  // if (m_TextAnnotations.size() != 3)
  // {
  //   m_TextAnnotations.clear();
  //   for (int i = 0; i < 3; i++)
  //   {
  //     m_TextAnnotations.push_back(mitk::TextAnnotation2D::New());
  //     auto renderer = GetRenderWindowPart()->GetQmitkRenderWindow(windownames[i].c_str())->GetRenderer();
  //     mitk::LayoutAnnotationRenderer::AddAnnotation(m_TextAnnotations.back(), renderer);
  //     m_TextAnnotations.back()->SetFontSize(15);
  //     float color[] = {0.7, 0.7, 0.7};
  //     m_TextAnnotations.back()->SetFontSize(20);

  //     m_TextAnnotations.back()->SetColor(color);
  //   }
  // }
  // for (auto anno : m_TextAnnotations)
  // {
  //   anno->SetText(text);
  // }
}

mitk::DataNode::Pointer m2Data::FindChildNodeRegex(mitk::DataNode::Pointer &parent, std::string regexString)
{
  auto deriv = this->GetDataStorage()->GetDerivations(parent.GetPointer(), nullptr, true);
  try
  {
    for (auto p : *deriv)
    {
      if (std::regex_match(p->GetName(), std::regex{regexString.c_str()}))
        return p;
    }
  }
  catch (std::exception &e)
  {
    MITK_WARN << e.what();
  }
  return nullptr;
}

void m2Data::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*part*/, const QList<mitk::DataNode::Pointer> &nodes)
{
  if (!nodes.empty())
  {
    auto node = nodes.front();
    if (auto image = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
    {
      QString labelText = str(boost::format("%.2f +/- %.2f Da") % image->GetCurrentX() % image->GetTolerance()).c_str();

      labelText += "\n";
      if (nodes.size() == 1)
        labelText += node->GetName().c_str();

      this->UpdateTextAnnotations(labelText.toStdString());
    }
  }
}

void m2Data::UpdateLevelWindow(const mitk::DataNode *node)
{
  if (auto msImageBase = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
  {
    mitk::LevelWindow lw;
    node->GetLevelWindow(lw);
    lw.SetAuto(msImageBase);
    if (m_Controls.CBUseFixedLevel->isChecked())
    {
      lw.SetLevelWindow(m_Controls.spnBxLevel->value(), m_Controls.spnBxWindow->value());
    }
    const_cast<mitk::DataNode *>(node)->SetLevelWindow(lw);
  }
}

void m2Data::NodeAdded(const mitk::DataNode *node)
{
  if (dynamic_cast<m2::OpenSlideImageIOHelperObject *>(node->GetData()))
  {
    OpenSlideImageNodeAdded(node);
  }
  else if (auto data = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
  {
    // !! Primary initialization !!
    this->ApplySettingsToImage(data);
    data->InitializeImageAccess();
    // this->RequestRenderWindowUpdate();

    SpectrumImageNodeAdded(node);
    auto xs = data->GetXAxis();

    if(m_Controls.spnBxMz->value() == 0)
      this->OnGenerateImageData(xs[xs.size()/2], FROM_GUI);
    else
      this->OnGenerateImageData(m_Controls.spnBxMz->value(), FROM_GUI);
  }
}

void m2Data::OpenSlideImageNodeAdded(const mitk::DataNode *node)
{
  if (auto openSlideIOHelper = dynamic_cast<m2::OpenSlideImageIOHelperObject *>(node->GetData()))
  {
    const auto name = node->GetName();
    auto dialog = new Qm2OpenSlideImageIOHelperDialog(m_Parent);
    dialog->SetOpenSlideImageIOHelperObject(openSlideIOHelper);
    auto result = dialog->exec();
    if (result == QDialog::Accepted)
    {
      try
      {
        auto data = dialog->GetData();
        // auto preview = dialog->GetPreviwData();
        // mitk::DataNode::Pointer parent = nullptr;
        // if (preview)
        // {
        // parent = mitk::DataNode::New();
        // parent->SetName(
        //   itksys::SystemTools::GetFilenameWithoutExtension(openSlideIOHelper->GetOpenSlideIO()->GetFileName()));
        // parent->SetData(dialog->GetPreviwData());
        // this->GetDataStorage()->Add(parent);
        // }
        if (data.size() == 1)
        {
          auto filter = m2::SubdivideImage2DFilter::New();
          filter->SetInput(data.back());
          filter->SetTileHeight((unsigned int)(1) << 13);
          filter->SetTileWidth((unsigned int)(1) << 13);
          filter->Update();

          const auto nX = filter->GetNumberOfTilesInX();
          const auto nY = filter->GetNumberOfTilesInY();
          MITK_INFO << nX << " " << nY;

          unsigned int k = 0;
          for (auto I : filter->GetOutputs())
          {
            auto node = mitk::DataNode::New();
            node->SetData(I);
            node->SetName(name + "_tile_" + std::to_string(k));
            node->SetVisibility(k < 2, nullptr);
            this->GetDataStorage()->Add(node);
            ++k;
          }
        }
        else
        {
          unsigned int i = 0;
          for (auto &I : data)
          {
            auto node = mitk::DataNode::New();
            node->SetData(I);
            node->SetName(
              itksys::SystemTools::GetFilenameWithoutExtension(openSlideIOHelper->GetOpenSlideIO()->GetFileName()) +
              "_" + std::to_string(i++));
            // node->SetVisibility(false, nullptr);
            this->GetDataStorage()->Add(node);
          }
        }
      }
      catch (std::exception &e)
      {
        MITK_ERROR << "Rendering Error" << e.what();
      }
    }
    // remove IO helper object from DS
    GetDataStorage()->Remove(node);
    delete dialog;
  }
}

void m2Data::ImzMLImageNodeAdded(const mitk::DataNode *) {}

void m2Data::FsmImageNodeAdded(const mitk::DataNode *)
{
  //
}

void m2Data::SpectrumImageNodeAdded(const mitk::DataNode *node)
{


  auto nodes = GetDataStorage()->GetAll();
  // resolve name-conflicts
  if (std::any_of(nodes->begin(),
                  nodes->end(),
                  [node](auto f) { return f != node && f->GetName().compare(node->GetName()) == 0; }))
  {
    bool ok;
    QString text = QInputDialog::getText(
      this->m_Parent, tr("Name conflict"), tr("Postfix:"), QLineEdit::Normal, QDir::home().dirName(), &ok);
    const_cast<mitk::DataNode *>(node)->SetName((node->GetName() + "_" + text.toStdString()).c_str());
  }

  auto lut = mitk::LookupTable::New();
  lut->SetType(mitk::LookupTable::LookupTableType::CIVIDS_TRANSPARENT);
  const_cast<mitk::DataNode *>(node)->SetProperty("LookupTable", mitk::LookupTableProperty::New(lut));
  const_cast<mitk::DataNode *>(node)->SetBoolProperty("binary", false);

  // add nodes to data storage (default: helper objects)
  if (auto spectrumImage = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
  {
    // -------------- add data interactor --------------

    auto interactor = m2::SpectrumImageDataInteractor::New();
    interactor->LoadStateMachine("PointSet.xml");
    interactor->SetEventConfig("PointSetConfig.xml");
    interactor->EnableInteraction(false);
    interactor->SetDataNode(const_cast<mitk::DataNode *>(node));
    interactor->SetDataStorage(GetDataStorage());
    interactor->EnableInteraction(true);
    const_cast<mitk::DataNode *>(node)->SetDataInteractor(interactor);

    m2::DefaultNodeProperties(node);

    // -------------- add Mask to datastorage --------------
    auto helperNode = mitk::DataNode::New();
    helperNode->SetName("MaskImage");
    helperNode->SetVisibility(m_Controls.showMaskImages->isChecked());
    helperNode->SetData(spectrumImage->GetMaskImage());
    helperNode->SetStringProperty("m2aia.helper.image.name", "MaskImage");

    // add hidden to DS
    helperNode->SetBoolProperty("helper object", true);
    this->GetDataStorage()->Add(helperNode, const_cast<mitk::DataNode *>(node));
    // consideration of the check boxes
    emit m_Controls.showMaskImages->toggled(m_Controls.showMaskImages->isChecked());

    // -------------- add ShiftImage to datastorage --------------
    if(spectrumImage->GetShiftImage()){

      helperNode = mitk::DataNode::New();
      helperNode->SetName("ShiftImage");
      helperNode->SetVisibility(false);
      helperNode->SetData(spectrumImage->GetShiftImage());
      helperNode->SetStringProperty("m2aia.helper.image.name", "ShiftImage");
      this->GetDataStorage()->Add(helperNode, const_cast<mitk::DataNode *>(node));
    }

    std::string inputLocation;
    node->GetStringProperty("MITK.IO.reader.inputlocation", inputLocation);
    auto fileName = itksys::SystemTools::GetFilenamePath(inputLocation) + "/" +
    itksys::SystemTools::GetFilenameWithoutLastExtension(inputLocation) + ".tSNE.nrrd";
    if(itksys::SystemTools::FileExists(fileName)){
      auto tsneImages = mitk::IOUtil::Load(fileName);
      helperNode = mitk::DataNode::New();
      helperNode->SetName("tSNE");
      helperNode->SetVisibility(false);
      helperNode->SetData(tsneImages[0]);
      helperNode->SetStringProperty("m2aia.helper.image.name", "tSNEImage");
      this->GetDataStorage()->Add(helperNode, const_cast<mitk::DataNode *>(node));
    }
    
    fileName = itksys::SystemTools::GetFilenamePath(inputLocation) + "/" +
    itksys::SystemTools::GetFilenameWithoutLastExtension(inputLocation) + ".PCA.nrrd";
    if(itksys::SystemTools::FileExists(fileName)){
      auto pcaImages = mitk::IOUtil::Load(fileName);
      helperNode = mitk::DataNode::New();
      helperNode->SetName("PCA");
      helperNode->SetVisibility(false);
      helperNode->SetData(pcaImages[0]);
      helperNode->SetStringProperty("m2aia.helper.image.name", "PCAImage");
      this->GetDataStorage()->Add(helperNode, const_cast<mitk::DataNode *>(node));
    }


    // -------------- add Index to datastorage --------------
    helperNode = mitk::DataNode::New();
    helperNode->SetName("IndexImage");
    helperNode->SetVisibility(m_Controls.showIndexImages->isChecked());
    helperNode->SetData(spectrumImage->GetIndexImage());
    helperNode->SetStringProperty("m2aia.helper.image.name", "IndexImage");
    helperNode->SetBoolProperty("binary", false);
    // add hidden to DS
    helperNode->SetBoolProperty("helper object", true);
    this->GetDataStorage()->Add(helperNode, const_cast<mitk::DataNode *>(node));
    // consideration of the check boxes
    emit m_Controls.showIndexImages->toggled(m_Controls.showIndexImages->isChecked());

    // -------------- add Normalization to datastorage --------------
    for (auto type : m2::NormalizationStrategyTypeList)
    {

      helperNode = mitk::DataNode::New();
      helperNode->SetName("NormalizationImage" + m2::to_string(type));
      helperNode->SetBoolProperty("binary", false);
      
      // here we add only the template images
      // initialization is done by UI interaction
      
      auto typeName = m2::NormalizationStrategyTypeNames[to_underlying(type)];
      fileName = itksys::SystemTools::GetFilenamePath(inputLocation) + "/" +
      itksys::SystemTools::GetFilenameWithoutLastExtension(inputLocation) + "." + typeName + ".nrrd";
      
      if(itksys::SystemTools::FileExists(fileName)){ 
        auto dataVector = mitk::IOUtil::Load(fileName);
        auto externalImage = dynamic_cast<mitk::Image *>(dataVector[0].GetPointer());
        spectrumImage->SetNormalizationImage(externalImage, type);
        spectrumImage->SetNormalizationImageStatus(type, true);
      }
      
      auto image = spectrumImage->GetNormalizationImage(type);
      helperNode->SetData(image); 

      // add hidden to DS
      helperNode->SetBoolProperty("helper object", true);
      this->GetDataStorage()->Add(helperNode, const_cast<mitk::DataNode *>(node));
      auto name = m2::NormalizationStrategyTypeNames.at(to_underlying(type))+ "Image";
      helperNode->SetStringProperty("m2aia.helper.image.normalization.name", name.c_str());
      helperNode->SetIntProperty("m2aia.helper.image.normalization.type", to_underlying(type));

      // clear the image data if not already initialized
      if(image && !spectrumImage->GetNormalizationImageStatus(type)){
        auto imageSize = image->GetDimensions();
        mitk::ImagePixelWriteAccessor<m2::NormImagePixelType, 3> acc(image);
        std::memset(acc.GetData(), 0, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(m2::NormImagePixelType));
      }

      // consideration of the check boxes
      auto checkBox = m_Controls.settings->findChild<QCheckBox *>(("ckBoxNormalizationImage" + m2::to_string(type)).c_str());
      emit checkBox->toggled(checkBox->isChecked());
    }

    // -------------- add Spectra to datastorage --------------
    // color for plots in spectrum view

    const auto AddSpectrum = [&node, this](std::string name,
                                           m2::SpectrumFormat type,
                                           std::string info,
                                           const std::vector<double> xs,
                                           const std::vector<double> ys,
                                           bool checkState,
                                           float alpha = 1.0)
    {
      auto intervalsNode = mitk::DataNode::New();
      auto intervals = m2::IntervalVector::New();
      intervals->SetType(type);
      intervals->SetInfo(info);
      intervals->SetProperty("m2aia.helper.spectrum.xaxis.count", mitk::IntProperty::New(xs.size()));

      if (auto image = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
        intervals->SetProperty("m2aia.image.pixel.count", mitk::IntProperty::New(image->GetNumberOfValidPixels()));

      using namespace std;
      auto &i = intervals->GetIntervals();
      transform(begin(xs), end(xs), begin(ys), back_inserter(i), [](auto &a, auto &b) { return m2::Interval(a, b); });
      
      intervalsNode->SetData(intervals);
      intervalsNode->SetName(name);
      intervalsNode->SetVisibility(this->Controls()->visibleOnInitialization->isChecked() && !m_ResetPreventDataStorageOverload.isRunning());
      intervalsNode->SetBoolProperty("helper object", !checkState);
      m2::CopyNodeProperties(node, intervalsNode);
      if ((unsigned int)(type) & (unsigned int)(m2::SpectrumFormat::Centroid))
        intervalsNode->SetOpacity(1 - alpha);
      intervalsNode->Modified();
      this->GetDataStorage()->Add(intervalsNode, const_cast<mitk::DataNode *>(node));     
    
    };

    const auto &xs = spectrumImage->GetXAxis();

    if (spectrumImage->GetSpectrumType().Format == m2::SpectrumFormat::ContinuousProfile ||
        spectrumImage->GetSpectrumType().Format == m2::SpectrumFormat::ProcessedProfile)
    {
      AddSpectrum("MaxSpectrum (" + node->GetName() + ")",
                  m2::SpectrumFormat::Profile,
                  "overview.max",
                  xs,
                  spectrumImage->GetSkylineSpectrum(),
                  m_Controls.showMaxSpectrum->isChecked());
      AddSpectrum("MeanSpectrum (" + node->GetName() + ")",
                  m2::SpectrumFormat::Profile,
                  "overview.mean",
                  xs,
                  spectrumImage->GetMeanSpectrum(),
                  m_Controls.showMeanSpectrum->isChecked());
      AddSpectrum("SingleSpectrum (" + node->GetName() + ")",
                  m2::SpectrumFormat::Profile,
                  "overview.single",
                  {},
                  {},
                  m_Controls.showSingleSpectrum->isChecked());
    }
    else
    {
      AddSpectrum("CentroidSpectrum (" + node->GetName() + ")",
                  m2::SpectrumFormat::Centroid,
                  "overview.centroids",
                  xs,
                  spectrumImage->GetMeanSpectrum(),
                  m_Controls.showCentroidSpectrum->isChecked(),
                  0.5);
      AddSpectrum("SingleSpectrum (" + node->GetName() + ")",
                  m2::SpectrumFormat::Centroid,
                  "overview.single",
                  {},
                  {},
                  m_Controls.showSingleSpectrum->isChecked(),
                  0.5);
    }

    m_ResetPreventDataStorageOverload.cancel();
    m_ResetPreventDataStorageOverload.setFuture(QtConcurrent::run([](){
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));}));
  }
}

void m2Data::NodeRemoved(const mitk::DataNode *node)
{
  if (dynamic_cast<m2::SpectrumImage *>(node->GetData()))
  {
    auto derivations = this->GetDataStorage()->GetDerivations(node);
    for (auto &&d : *derivations)
    {
      this->GetDataStorage()->Remove(d);
    }
  }
}
