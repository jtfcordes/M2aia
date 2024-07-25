/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/jtfcordes/m2aia for details.

===================================================================*/

#ifndef __M2_BROWSER_PREFERENCE_PAGE_H
#define __M2_BROWSER_PREFERENCE_PAGE_H

#include "berryIQtPreferencePage.h"
#include <QProcess>
#include <QString>
#include <mitkIPreferences.h>

class QWidget;
class QCheckBox;
class QComboBox;
class QSpinBox;
class ctkDirectoryButton;

class Ui_m2PreferencePage;

/**
 * \class m2BrowserPreferencesPage
 * \brief Preference page for the MatchPoint Browser plugin
 */
class m2BrowserPreferencesPage : public QObject, public berry::IQtPreferencePage
{
    Q_OBJECT
    Q_INTERFACES(berry::IPreferencePage)

public:
    m2BrowserPreferencesPage();
    ~m2BrowserPreferencesPage() override;

private:
	QProcess *m_ElastixProcess;
	std::string m_ElastixPath;

	QProcess *m_TransformixProcess;
	std::string m_TransformixPath;

public slots:
	void OnElastixButtonClicked();
	void OnElastixProcessError(QProcess::ProcessError error);
	void OnElastixProcessFinished(int exitCode, QProcess::ExitStatus exitStatus);

	void OnTransformixButtonClicked();
	void OnTransformixProcessError(QProcess::ProcessError error);
	void OnTransformixProcessFinished(int exitCode, QProcess::ExitStatus exitStatus);

	void OnBinsSpinBoxValueChanged(int);

	void OnUseMinIntensity(bool v);
	void OnUseMaxIntensity(bool v);
	void OnUseSamplingPoints(bool v);
	void OnUseMinimalImagingArea(bool v);

	void CreateQtControl(QWidget* parent) override;
	QWidget* GetQtControl() const override;
	void Init(berry::IWorkbench::Pointer) override;
	void PerformCancel() override;
	bool PerformOk() override;
	void Update() override;

private:

    QString ConvertToString(const QStringList& list);
	mitk::IPreferences* m_Preferences;
	QScopedPointer<Ui_m2PreferencePage> m_Ui;
	QWidget* m_Control;
};

#endif
