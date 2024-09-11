/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/
#pragma once

#include <QmitkAbstractDataNodeAction.h>
#include <QAction>
#include <mitkImage.h>


class QmitkDataNodeExportComponentAction : public QAction, public QmitkAbstractDataNodeAction
{
  Q_OBJECT

public:

  QmitkDataNodeExportComponentAction(QWidget* parent, berry::IWorkbenchPartSite::Pointer workbenchPartSite);
  QmitkDataNodeExportComponentAction(QWidget* parent = nullptr, berry::IWorkbenchPartSite* workbenchPartSite = nullptr);

  static mitk::Image::Pointer ExportComponentImage(mitk::Image * img, unsigned int i);

private Q_SLOTS:

  void OnActionTriggered(bool);

protected:

  void InitializeAction() override;

};
