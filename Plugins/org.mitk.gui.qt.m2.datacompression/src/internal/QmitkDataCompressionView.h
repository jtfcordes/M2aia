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

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include <QmitkSingleNodeSelectionWidget.h>
#include <ui_DataCompressionViewControls.h>
#include <mitkImage.h>


class QmitkDataCompressionView : public QmitkAbstractView
{
  Q_OBJECT

public:
  static const std::string VIEW_ID;
  void CreateQtPartControl(QWidget* parent) override;

private slots:
  void OnStartPCA();
  void OnStartTSNE();
  void OnStartKMeans();
  void OnSaveDataCompressionResults();
  

private:
  mitk::Image::Pointer ResampleVectorImage(mitk::Image::Pointer lowResImage, mitk::Image::Pointer referenceImage);
  void SetFocus() override;
  QWidget * m_Parent;

  Ui::DataCompressionViewControls m_Controls;
};

