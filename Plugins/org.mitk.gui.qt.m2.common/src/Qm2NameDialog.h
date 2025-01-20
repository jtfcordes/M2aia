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
#include <org_mitk_gui_qt_m2_common_Export.h>
#include "ui_Qm2NameDialogControls.h"
#include <QDialog>

class MITK_M2_CORE_HELPER_EXPORT Qm2NameDialog : public QDialog
{
  Q_OBJECT

public:
  explicit Qm2NameDialog(QWidget *parent = nullptr) : QDialog(parent)
  {
    m_Controls.setupUi(this);
    this->setSizeGripEnabled(true);
  }

  virtual  ~Qm2NameDialog(){}

  int exec() override
  {
    return QDialog::exec();
  }

  std::string GetName() { return m_Controls.lineEdit->text().toStdString(); }
  void SetName(std::string v) { return m_Controls.lineEdit->setText(v.c_str()); }
  
private:
  
  Ui::NameDialog m_Controls;
 
  
};
