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
#include <QThreadPool>
#include <QmitkAbstractView.h>
#include <QmitkDataStorageComboBox.h>
#include <QmitkDataStorageComboBoxWithSelectNone.h>
#include <QmitkPointListWidget.h>
#include <berryISelectionListener.h>
#include <itkVectorContainer.h>
#include <m2CommunicationService.h>
#include <m2ImzMLMassSpecImage.h>
#include <m2MSImageBase.h>
#include <mitkColorBarAnnotation.h>
#include <mitkColorSequenceRainbow.h>
#include <mitkDataNode.h>
#include <mitkNodePredicateDataType.h>
#include <mitkTextAnnotation2D.h>

/**
  \brief ImsPosition

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
namespace mitk
{
  class ImzMLMassSpecImage;
}

class QmitkMultiNodeSelectionWidget;

class m2Data : public QmitkAbstractView
{
  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

public:
  static const std::string VIEW_ID;
  using NodesVectorType = m2::CommunicationService::NodesVectorType;

  /**
   * @brief Get the Overview Spectrum Type object
   *
   * @return MSImageBase::OverviewSpectrumType
   */
  m2::OverviewSpectrumType GetOverviewSpectrumType() { return m_CurrentOverviewSpectrumType; }

  /**
   * @brief
   *
   * @return Ui::imsDataControls*
   */
  Ui::imsDataControls *Controls() { return &m_Controls; }


  /**
   * @brief Return all nodes selected by GUI-defined condistions.
   *
   * @return A set of data nodes.
   *
   */
  NodesVectorType::Pointer AllNodes();

  /**
   * @brief Apply settings to nodes.
   *
   * @param v A set of data node of mass spec images.
   * Configures the I/O strategy state of the data objects.
   */
  void ApplySettingsToNodes(NodesVectorType::Pointer v);

  /**
   * @brief Apply settings to a single MS image.
   * Configures the I/O strategy state of a single MS image.
   * @param image A mass spectrometry image.
   */
  void ApplySettingsToImage(m2::ImzMLMassSpecImage *image);

  m2::NormalizationStrategyType GuiToNormalizationStrategyType();
  m2::IonImageGrabStrategyType GuiToIonImageGrabStrategyType();

  m2::SmoothingType GuiToSmoothingStrategyType();

  m2::BaselineCorrectionType GuiToBaselineCorrectionStrategyType();

public slots:
  void UpdateColorBarAndRenderWindows();
  void OnEqualizeLW();
  void OnApplyTiling();
  void OnResetTiling();
  void OnNextIonImage();
  void OnPrevIonImage();
  void OnProcessingNodesRequested(const QString &);
  void EmitIonImageReference();

  /**
   * @brief GUI side invocation of an ion image grab.
   * Update render window text annotation.
   * Configures I/O strategy state of the mass spec image.
   *
   * @param mz Center of the mass range window
   * @param tol Tolerance (tol) to control mass range window that is concidred.
   */
  void OnGrabIonImage(qreal mz, qreal tol);
  void OnGrabIonImageFinished(mitk::DataNode * parent, mitk::Image * image);

signals:
  void SendNodes(NodesVectorType::Pointer vec);

protected:

	mitk::DataNode::Pointer FindChildNodeWithNameContaining(mitk::DataNode::Pointer & parent, const std::string & subStr);


  virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer part,
                                  const QList<mitk::DataNode::Pointer> &nodes) override;
  void UpdateLevelWindow(const mitk::DataNode *node);
  virtual void NodeAdded(const mitk::DataNode *node) override;
  virtual void NodeRemoved(const mitk::DataNode *node) override;
  virtual void CreateQtPartControl(QWidget *parent) override;
  virtual void SetFocus() override {}

  // 20201023: custom selection service did not work as expected
  //  m2::SelectionProvider::Pointer m_SelectionProvider;
  // void SetSelectionProvider() override;

  // 20201023: custom selection service did not work as expected
  // QScopedPointer<berry::ISelectionListener> m_CustomSelectionChangedListener, m_NodeSelectionListener;
  /*virtual void OnCustomSelectionChanged(const berry::IWorkbenchPart::Pointer &sourcepart,
                                        const berry::ISelection::ConstPointer &selection);*/

  Ui::imsDataControls m_Controls;
  QThreadPool m_pool;
  m2::OverviewSpectrumType m_CurrentOverviewSpectrumType = m2::OverviewSpectrumType::Maximum;

  m2::IonImageReference::Pointer m_IonImageReference;

  mitk::ColorSequenceRainbow m_ColorSequence;

  /*!
   * Main element holding a list of DataNodes containing MassSpecBaseData objects.
   * Thes selected node in the combobox may be in focus of the processing controlled by this view.
   */
  QmitkMultiNodeSelectionWidget *m_MassSpecDataNodeSelectionWidget;

  /*!
   * Main element holding a list of DataNodes containing MassSpecBaseData objects.
   * Thes selected node in the combobox may be in focus of the processing controlled by this view.
   */
  QmitkMultiNodeSelectionWidget *m_PointSetDataNodeSelectionWidget;

  QmitkPointListWidget *m_PointListWidget;

  std::vector<mitk::TextAnnotation2D::Pointer> m_TextAnnotations;
  std::vector<mitk::ColorBarAnnotation::Pointer> m_ColorBarAnnotations;
  void UpdateTextAnnotations(std::string text);
};
