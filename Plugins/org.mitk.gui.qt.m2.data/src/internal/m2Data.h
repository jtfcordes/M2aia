/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/jtfcordes/m2aia for details.

===================================================================*/

#pragma once
#include "ui_m2Data.h"
#include <QmitkAbstractView.h>

#include <mitkImage.h>
#include <mitkTextAnnotation2D.h>
#include <mitkColorBarAnnotation.h>

#include <m2UIUtils.h>

// #include <QThreadPool>
#include <QtConcurrent>
/**
  \brief DataView

  \warning  The main view for MSI data processing. 

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
namespace m2
{
  class SpectrumImage;
  class ImzMLSpectrumImage;
}

class QmitkMultiNodeSelectionWidget;

class m2Data : public QmitkAbstractView
{
  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

public:
  static const std::string VIEW_ID;
  

  /**
   * @brief Get the Overview Spectrum Type object
   *
   * @return SpectrumImage::OverviewSpectrumType
   */
  m2::SpectrumType GetOverviewSpectrumType() { return m_CurrentOverviewSpectrumType; }

  /**
   * @brief
   *
   * @return Ui::imsDataControls*
   */
  Ui::imsDataControls *Controls() { return &m_Controls; }


  /**
   * @brief Apply settings to nodes.
   *
   * @param v A set of data node of mass spec images.
   * Configures the I/O strategy state of the data objects.
   */
  void ApplySettingsToNodes(m2::UIUtils::NodesVectorType::Pointer v);

  /**
   * @brief Apply settings to a single MS image.
   * Configures the I/O strategy state of a single MS image.
   * @param image A mass spectrometry image.
   */
  void ApplySettingsToImage(m2::SpectrumImage *image);

  void InitNormalizationControls();
  void InitIntensityTransformationControls();
  void InitRangePoolingControls();
  void InitSmoothingControls();
  void InitToleranceControls();
  void InitBaselineCorrectionControls();
  void InitImageNormalizationControls();
  void InitImageSmoothingControls();
  
  

public slots:
  void OnCreateNextImage();
  void OnCreatePrevImage();
  void OnCreatePrevPeakImage();
  void OnCreateNextPeakImage();
  void OnCreateShiftMap();
  void OnIncreaseTolerance();
  void OnDecreaseTolerance();
  void OnRenderSpectrumImages(double min, double max);
  void OnGenerateImageData(mitk::DataNode::Pointer node, qreal mz, qreal tol, bool emitRangeChanged);
  void OnGenerateImageData(qreal mz, qreal tol);
  
  m2::NormalizationStrategyType GuiToNormalizationStrategyType();
  m2::IntensityTransformationType GuiToIntensityTransformationStrategyType();
  m2::RangePoolingStrategyType GuiToRangePoolingStrategyType();  
  m2::SmoothingType GuiToSmoothingStrategyType();  
  m2::BaselineCorrectionType GuiToBaselineCorrectionStrategyType();
  m2::ImageNormalizationStrategyType GuiToImageNormalizationStrategyType();
  m2::ImageSmoothingStrategyType GuiToImageSmoothingStrategyType();
  

signals:
  void SendNodes(m2::UIUtils::NodesVectorType::Pointer vec);

protected:
  void CreateQtPartControl(QWidget *parent) override;
  void NodeAdded(const mitk::DataNode *node) override;
  QFutureWatcher<void> m_ResetPreventDataStorageOverload;
  void NodeRemoved(const mitk::DataNode *node) override;
  void SetFocus() override {}

	mitk::DataNode::Pointer FindChildNodeRegex(mitk::DataNode::Pointer & parent, std::string regexString);

  virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer part,
                                  const QList<mitk::DataNode::Pointer> &nodes) override;
  void UpdateLevelWindow(const mitk::DataNode *node);
  void UpdateSpectrumImageTable(const mitk::DataNode* node);
  void OpenSlideImageNodeAdded(const mitk::DataNode *node);
  void ImzMLImageNodeAdded(const mitk::DataNode *node);
  void FsmImageNodeAdded(const mitk::DataNode *node);
  void SpectrumImageNodeAdded(const mitk::DataNode *node);

  // 20201023: custom selection service did not work as expected
  //  m2::SelectionProvider::Pointer m_SelectionProvider;
  // void SetSelectionProvider() override;

  // 20201023: custom selection service did not work as expected
  // QScopedPointer<berry::ISelectionListener> m_CustomSelectionChangedListener, m_NodeSelectionListener;
  /*virtual void OnCustomSelectionChanged(const berry::IWorkbenchPart::Pointer &sourcepart,
                                        const berry::ISelection::ConstPointer &selection);*/

  Ui::imsDataControls m_Controls;
  QWidget * m_Parent = nullptr;
  bool m_InitializeNewNode = false;

  QThreadPool m_pool;
  m2::SpectrumType m_CurrentOverviewSpectrumType = m2::SpectrumType::Maximum;


  // m2::IonImageReference::Pointer m_IonImageReference;
  /*!
   * Main element holding a list of DataNodes containing PointSet objects.
   * Thes selected node in the combobox may be in focus of the processing controlled by this view.
   */
  QmitkMultiNodeSelectionWidget *m_PointSetDataNodeSelectionWidget;

  // QmitkPointListWidget *m_PointListWidget;

  std::vector<mitk::TextAnnotation2D::Pointer> m_TextAnnotations;
  std::vector<mitk::ColorBarAnnotation::Pointer> m_ColorBarAnnotations;
  void UpdateTextAnnotations(std::string text);

  const int FROM_GUI = -1;

  mitk::IPreferences * m_M2aiaPreferences;
};
