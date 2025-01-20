/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/jtfcordes/m2aia for details.

===================================================================*/

#ifndef m2Reconstruction3D_h
#define m2Reconstruction3D_h

#include "ui_m2Reconstruction3DControls.h"
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QFuture>
#include <QMimeData>
#include <QPoint>
#include <QPointF>
#include <QmitkAbstractView.h>
#include <berryISelectionListener.h>
#include <m2ImzMLSpectrumImage.h>
#include <m2SpectrumImage.h>
#include <m2SpectrumImageStack.h>
#include <mitkImage.h>
#include <mitkPointSet.h>
#include <qfuturewatcher.h>
#include <qlistwidget.h>
#include <qprocess.h>
#include <ui_ParameterFileEditorDialog.h>

/**
  \brief m2Reconstruction3D

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class m2Reconstruction3D : public QmitkAbstractView
{
  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

public:
  static const std::string VIEW_ID;
  QFutureWatcher<void> m_ExportProcessFutureWatcher;
  QFutureWatcher<void> m_ReconstructionFutureWatcher;

protected:
  virtual void CreateQtPartControl(QWidget *parent) override;
  virtual void SetFocus() override {}
  QWidget *m_Parent;

  struct DataTuple
  {
    mitk::Image::Pointer image, mask;
    mitk::DataNode::Pointer node;
    mitk::PointSet::Pointer points;
  };

  Ui::m2Reconstruction3DControls m_Controls;
  std::shared_ptr<QProcess> m_Process;

  DataTuple GetImageDataById(unsigned int id, QListWidget *listWidget);

  std::shared_ptr<m2::ElxRegistrationHelper> RegistrationStep(unsigned int fixedId,
                                                                      QListWidget *fixedSource,
                                                                      std::shared_ptr<m2::ElxRegistrationHelper> fixedTransformer,
                                                                      unsigned int movingId,
                                                                      QListWidget *movingSource);
  std::vector<std::string> GetParameters();

  QListWidget *m_List1, *m_List2;

  std::map<unsigned int, DataTuple> m_referenceMap;

  Ui::elxParameterFileEditor m_ParameterFileEditorControls;
  QDialog *m_ParameterFileEditor;
  std::vector<std::string> m_ParameterFiles;

public slots:
  void OnUpdateList();
  void OnStartStacking();
  void OnStartExport();
};

#endif // m2Reconstruction3D_h
