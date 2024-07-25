/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/jtfcordes/m2aia for details.

===================================================================*/

#include "m2BrowserPreferencesPage.h"

#include <mitkCoreServices.h>
#include <mitkIPreferencesService.h>
#include <mitkIPreferences.h>

#include <mitkExceptionMacro.h>
#include <ui_m2PreferencePage.h>

#include <QFileDialog>
#include <QTextCodec>

static mitk::IPreferences* GetPreferences()
{
	auto* preferencesService = mitk::CoreServices::GetPreferencesService();
	return preferencesService->GetSystemPreferences();
}

m2BrowserPreferencesPage::m2BrowserPreferencesPage()
    :
	m_ElastixProcess(nullptr),
    m_TransformixProcess(nullptr),
    m_Preferences(GetPreferences()),
    m_Ui(new Ui_m2PreferencePage),
    m_Control(nullptr)

{
}

m2BrowserPreferencesPage::~m2BrowserPreferencesPage()
{
}

void m2BrowserPreferencesPage::CreateQtControl(QWidget* parent)
{
	m_Control = new QWidget(parent);
	m_ElastixProcess = new QProcess(m_Control);
	m_TransformixProcess = new QProcess(m_Control);

	m_Ui->setupUi(m_Control);

	// init widgets	
	auto bins = m_Preferences->GetInt("m2aia.view.spectrum.bins", 15000);
	m_Ui->spnBxBins->setValue(bins);

	m_Ui->useMaxIntensity->setChecked(m_Preferences->GetBool("m2aia.view.spectrum.useMaxIntensity", true));
    m_Ui->useMinIntensity->setChecked(m_Preferences->GetBool("m2aia.view.spectrum.useMinIntensity", true));
	
	// auto v = ;
	// m_Preferences->PutBool("m2aia.view.spectrum.showSamplingPoints",v);
	m_Ui->showSamplingPoints->setChecked(m_Preferences->GetBool("m2aia.view.spectrum.showSamplingPoints", false));
	m_Ui->minimalImagingArea->setChecked(m_Preferences->GetBool("m2aia.view.image.minimal_area", true));

	m_ElastixPath = m_Preferences->Get("m2aia.external.elastix", "");
	if (!m_ElastixPath.empty())
		m_ElastixProcess->start(m_ElastixPath.c_str(), QStringList() << "--version", QProcess::ReadOnly);

	m_TransformixPath = m_Preferences->Get("m2aia.external.transformix", "");
	if (!m_TransformixPath.empty())
		m_TransformixProcess->start(m_TransformixPath.c_str(), QStringList() << "--version", QProcess::ReadOnly);

	// connect widgts
	connect(m_ElastixProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(OnElastixProcessError(QProcess::ProcessError)));
	connect(m_ElastixProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(OnElastixProcessFinished(int, QProcess::ExitStatus)));
	connect(m_Ui->elastixButton, SIGNAL(clicked()), this, SLOT(OnElastixButtonClicked()));

	connect(
		m_TransformixProcess, SIGNAL(error(QProcess::ProcessError)), this, SLOT(OnTransformixProcessError(QProcess::ProcessError)));
	connect(m_TransformixProcess,
		SIGNAL(finished(int, QProcess::ExitStatus)),
		this,
		SLOT(OnTransformixProcessFinished(int, QProcess::ExitStatus)));
	connect(m_Ui->transformixButton, SIGNAL(clicked()), this, SLOT(OnTransformixButtonClicked()));

	connect(m_Ui->spnBxBins, SIGNAL(valueChanged(int)), this, SLOT(OnBinsSpinBoxValueChanged(int)));
	connect(m_Ui->useMaxIntensity, SIGNAL(toggled(bool)), this, SLOT(OnUseMaxIntensity(bool)));
	connect(m_Ui->useMinIntensity, SIGNAL(toggled(bool)), this, SLOT(OnUseMinIntensity(bool)));
	connect(m_Ui->minimalImagingArea, SIGNAL(toggled(bool)), this, SLOT(OnUseMinimalImagingArea(bool)));
	connect(m_Ui->showSamplingPoints, SIGNAL(toggled(bool)), this, SLOT(OnUseSamplingPoints(bool)));

}


void m2BrowserPreferencesPage::OnBinsSpinBoxValueChanged(int value){
	m_Preferences->PutInt("m2aia.view.spectrum.bins", value);
}


void m2BrowserPreferencesPage::OnElastixButtonClicked()
{
	QString filter = "elastix executable ";

#if defined(WIN32)
	filter += "(elastix.exe)";
#else
	filter += "(elastix)";
#endif

	QString path = QFileDialog::getOpenFileName(m_Control, "Elastix", "", filter);

	if (!path.isEmpty())
	{
		m_ElastixPath = path.toStdString();
		m_ElastixProcess->start(path, QStringList() << "--version", QProcess::ReadOnly);
	}
}

void m2BrowserPreferencesPage::OnElastixProcessError(QProcess::ProcessError)
{
	m_ElastixPath.clear();
	m_Ui->elastixLineEdit->clear();
}

void m2BrowserPreferencesPage::OnElastixProcessFinished(int exitCode, QProcess::ExitStatus exitStatus)
{
	if (exitStatus == QProcess::NormalExit && exitCode == 0)
	{
		QString output = QTextCodec::codecForName("UTF-8")->toUnicode(m_ElastixProcess->readAllStandardOutput());

		if (output.startsWith("elastix"))
		{
			m_Ui->elastixLineEdit->setText(m_ElastixPath.c_str());
			return;
		}
	}

	m_ElastixPath.clear();
	m_Ui->elastixLineEdit->clear();
}

void m2BrowserPreferencesPage::OnTransformixButtonClicked()
{
	QString filter = "transformix executable ";

#if defined(WIN32)
	filter += "(transformix.exe)";
#else
	filter += "(transformix)";
#endif

	QString path = QFileDialog::getOpenFileName(m_Control, "Transformix", "", filter);

	if (!path.isEmpty())
	{
		m_TransformixPath = path.toStdString();
		m_TransformixProcess->start(path, QStringList() << "--version", QProcess::ReadOnly);
	}
}

void m2BrowserPreferencesPage::OnTransformixProcessError(QProcess::ProcessError)
{
	m_TransformixPath.clear();
	m_Ui->transformixLineEdit->clear();
}

void m2BrowserPreferencesPage::OnTransformixProcessFinished(int exitCode, QProcess::ExitStatus exitStatus)
{
	if (exitStatus == QProcess::NormalExit && exitCode == 0)
	{
		QString output = QTextCodec::codecForName("UTF-8")->toUnicode(m_TransformixProcess->readAllStandardOutput());

		if (output.startsWith("transformix"))
		{
			m_Ui->transformixLineEdit->setText(m_TransformixPath.c_str());
			return;
		}
	}

	m_TransformixPath.clear();
	m_Ui->transformixLineEdit->clear();
}

QWidget* m2BrowserPreferencesPage::GetQtControl() const
{
	return m_Control;
}

void m2BrowserPreferencesPage::Init(berry::IWorkbench::Pointer)
{
	
}

void m2BrowserPreferencesPage::PerformCancel()
{
}

bool m2BrowserPreferencesPage::PerformOk()
{
	m_Preferences->Put("m2aia.external.elastix", m_ElastixPath);
	m_Preferences->Put("m2aia.external.transformix", m_TransformixPath);
	return true;
}

void m2BrowserPreferencesPage::OnUseSamplingPoints(bool v){
	m_Preferences->PutBool("m2aia.view.spectrum.showSamplingPoints", v);
}


void m2BrowserPreferencesPage::OnUseMaxIntensity(bool v){
	m_Preferences->PutBool("m2aia.view.spectrum.useMaxIntensity", v);
}

void m2BrowserPreferencesPage::OnUseMinIntensity(bool v){
	m_Preferences->PutBool("m2aia.view.spectrum.useMinIntensity", v);
}

void m2BrowserPreferencesPage::OnUseMinimalImagingArea(bool v){
	m_Preferences->PutBool("m2aia.view.image.minimal_area", v);
}

void m2BrowserPreferencesPage::Update()
{
	
	//optin
   
}
