/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/jtfcordes/m2aia for details.

===================================================================*/

#include "QmitkM2aiaAppWorkbenchAdvisor.h"

#include "internal/QmitkM2aiaApplicationPlugin.h"
#include <QApplication>
// #include <QDesktopWidget>
#include <QMainWindow>
#include <QMenuBar>
#include <QPoint>
#include <QRect>
#include <QmitkExtWorkbenchWindowAdvisor.h>
#include <QmitkM2aiaAboutDialog.h>
#include <QmitkM2aiaViewAction.h>
#include <mitkCoreServices.h>
#include <mitkIPreferences.h>
#include <mitkIPreferencesService.h>
#include <berryIWorkbenchWindow.h>
#include <berryPlatform.h>
#include <berryQtPreferences.h>
#include <berryWorkbenchPlugin.h>
#include <berryWorkbenchPreferenceConstants.h>
#include <mitkBaseRenderer.h>
#include <mitkNodePredicateData.h>
#include <mitkCameraController.h>
#include <mitkImage.h>
#include <mitkM2aiaVersion.h>
#include <mitkVersion.h>
#include <mitkIOUtil.h>
#include <mitkIDataStorageService.h>
#include <qaction.h>
#include <qevent.h>
#include <qtoolbar.h>
#include <qtoolbutton.h>

#include <internal/QmitkM2aiaApplicationPlugin.h>

QmitkM2aiaAppWorkbenchWindowAdvisor::QmitkM2aiaAppWorkbenchWindowAdvisor(
  berry::WorkbenchAdvisor *wbAdvisor, berry::IWorkbenchWindowConfigurer::Pointer configurer)
  : QmitkExtWorkbenchWindowAdvisor(wbAdvisor, configurer)
{
}

void QmitkM2aiaAppWorkbenchWindowAdvisor::loadDataFromDisk(const QStringList &arguments, bool globalReinit)
{
  if (!arguments.empty())
  {
    auto _context = QmitkM2aiaApplicationPlugin::GetDefault()->GetPluginContext();
    ctkServiceReference serviceRef = _context->getServiceReference<mitk::IDataStorageService>();
    if (serviceRef)
    {
      mitk::IDataStorageService* dataStorageService = _context->getService<mitk::IDataStorageService>(serviceRef);
      mitk::DataStorage::Pointer dataStorage = dataStorageService->GetDefaultDataStorage()->GetDataStorage();

      int argumentsAdded = 0;
      for (int i = 0; i < arguments.size(); ++i)
      {
        if (arguments[i].startsWith("--m2."))
        { //By convention no further files are specified as soon as a flow arguments comes.
          
          break;
        }
        else if(arguments[i].endsWith(".nrrd") || arguments[i].endsWith(".imzML", Qt::CaseSensitivity::CaseInsensitive))
        {
          try
          {
            std::string filename = arguments[i].toStdString();
            mitk::IOUtil::Load(filename);
            argumentsAdded++;
          }
          catch (...)
          {
            MITK_WARN << "Failed to load command line argument: " << arguments[i].toStdString();
          }
        }
      } // end for each command line argument

      if (argumentsAdded > 0 && globalReinit)
      {
        // calculate bounding geometry
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(dataStorage);
      }
    }
    else
    {
      MITK_ERROR << "A service reference for mitk::IDataStorageService does not exist";
    }
  }
}

void QmitkM2aiaAppWorkbenchWindowAdvisor::PostWindowOpen()
{
  QmitkExtWorkbenchWindowAdvisor::PostWindowOpen();
  berry::IWorkbenchWindowConfigurer::Pointer configurer = GetWindowConfigurer();
  configurer->GetWindow()->GetWorkbench()->GetIntroManager()->ShowIntro(configurer->GetWindow(), false);

  // very bad hack...
  berry::IWorkbenchWindow::Pointer window = this->GetWindowConfigurer()->GetWindow();
  QMainWindow *mainWindow = qobject_cast<QMainWindow *>(window->GetShell()->GetControl());
  mainWindow->setContextMenuPolicy(Qt::PreventContextMenu);

  
  // ==== RenderwindowViewNames ======================

  mitk::IPreferencesService *prefService = mitk::CoreServices::GetPreferencesService();
  Q_ASSERT(prefService);


  // ==== Application menu ============================
  QMenuBar *menuBar = mainWindow->menuBar();
  menuBar->setContextMenuPolicy(Qt::PreventContextMenu);
  auto actions = menuBar->actions();
  for (auto a : actions)
  {
    if (a->text() == "&Help")
    {
      auto helpMenu = a->menu();
      helpMenu->addAction("&About M²aia Workbench", this, []() { QmitkM2aiaAboutDialog().exec(); });
      for (auto b : helpMenu->actions())
      {
        if (b->text() == "&About")
        {
          b->setText("&About MITK");
        }
      }
    }
  }

  const std::vector<QString> viewCategories = {"Spectrum imaging", "Segmentation", "Spectrum imaging (Docker)"};

  {
    // auto prefService = berry::WorkbenchPlugin::GetDefault()->GetPreferencesService();
    // auto stylePrefs = prefService->GetSystemPreferences()->Node(berry::QtPreferences::QT_STYLES_NODE);
    // bool showCategoryNames = stylePrefs->GetBool(berry::QtPreferences::QT_SHOW_TOOLBAR_CATEGORY_NAMES, true);

    // Order view descriptors by category

    QMultiMap<QString, berry::IViewDescriptor::Pointer> categoryViewDescriptorMap;
    std::map<QString, berry::IViewDescriptor::Pointer> VDMap;
    // sort elements (converting vector to map...)
    QList<berry::IViewDescriptor::Pointer>::const_iterator iter;

    berry::IViewRegistry *viewRegistry = berry::PlatformUI::GetWorkbench()->GetViewRegistry();
    const QList<berry::IViewDescriptor::Pointer> viewDescriptors = viewRegistry->GetViews();

    // bool skip = false;
    for (iter = viewDescriptors.begin(); iter != viewDescriptors.end(); ++iter)
    {
      if ((*iter)->GetId() == "org.blueberry.ui.internal.introview")
        continue;
      if ((*iter)->GetId() == "org.mitk.views.imagenavigator")
        continue;
      if ((*iter)->GetId() == "org.mitk.views.viewnavigatorview")
        continue;
      if ((*iter)->GetId() == "org.mitk.views.m2.data")
        continue;
      if ((*iter)->GetId() == "org.mitk.views.m2.spectrum")
        continue;

      std::pair<QString, berry::IViewDescriptor::Pointer> p((*iter)->GetLabel(), (*iter));
      VDMap.insert(p);
    }
    for (auto labelViewDescriptorPair : VDMap)
    {
      auto viewDescriptor = labelViewDescriptorPair.second;
      auto category =
        !viewDescriptor->GetCategoryPath().isEmpty() ? viewDescriptor->GetCategoryPath().back() : QString();

      categoryViewDescriptorMap.insert(category, viewDescriptor);
    }

    // Create a separate toolbar for each category
    berry::IWorkbenchWindow::Pointer window = this->GetWindowConfigurer()->GetWindow();
    QMainWindow *mainWindow = qobject_cast<QMainWindow *>(window->GetShell()->GetControl());

    for (auto category : categoryViewDescriptorMap.uniqueKeys())
    {
      if (std::find(std::begin(viewCategories), std::end(viewCategories), category) == std::end(viewCategories))
        continue;

      auto viewDescriptorsInCurrentCategory = categoryViewDescriptorMap.values(category);

      if (!viewDescriptorsInCurrentCategory.isEmpty())
      {
        auto toolbar = new QToolBar;
        toolbar->setObjectName(category + " View Toolbar");
        mainWindow->addToolBar(toolbar);

        if (!category.isEmpty())
        {
          auto categoryButton = new QToolButton;
          categoryButton->setToolButtonStyle(Qt::ToolButtonTextOnly);
          categoryButton->setText("");
          categoryButton->setStyleSheet("background: transparent; margin: 0; padding: 0;");
          toolbar->addWidget(categoryButton);

          connect(
            categoryButton,
            &QToolButton::clicked,
            [toolbar]()
            {
              for (QWidget *widget : toolbar->findChildren<QWidget *>())
              {
                if (QStringLiteral("qt_toolbar_ext_button") == widget->objectName() && widget->isVisible())
                {
                  QMouseEvent pressEvent(
                    QEvent::MouseButtonPress, QPointF(0.0f, 0.0f), Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
                  QMouseEvent releaseEvent(
                    QEvent::MouseButtonRelease, QPointF(0.0f, 0.0f), Qt::LeftButton, Qt::LeftButton, Qt::NoModifier);
                  QApplication::sendEvent(widget, &pressEvent);
                  QApplication::sendEvent(widget, &releaseEvent);
                }
              }
            });
        }

        for (auto viewDescriptor : viewDescriptorsInCurrentCategory)
        {
          auto viewAction = new m2::QmitkM2aiaViewAction(window, viewDescriptor);
          toolbar->addAction(viewAction);
        }
      }
    }
  }

  {
    auto m_Preferences = prefService->GetSystemPreferences()->Node("org.mitk.editors.stdmultiwidget");
    m_Preferences->Put("DepartmentLogo", ":/M2aiaApplication/defaultWatermark.png");
  }


  loadDataFromDisk(berry::Platform::GetApplicationArgs(), true);
}

void QmitkM2aiaAppWorkbenchAdvisor::Initialize(berry::IWorkbenchConfigurer::Pointer configurer)
{
  berry::QtWorkbenchAdvisor::Initialize(configurer);

  configurer->SetSaveAndRestore(true);
}

berry::WorkbenchWindowAdvisor *QmitkM2aiaAppWorkbenchAdvisor::CreateWorkbenchWindowAdvisor(
  berry::IWorkbenchWindowConfigurer::Pointer configurer)
{
  QList<QString> perspExcludeList;
  perspExcludeList.push_back("org.blueberry.perspectives.help");

  QList<QString> viewExcludeList;
  viewExcludeList.push_back("org.mitk.views.controlvisualizationpropertiesview");
  viewExcludeList.push_back("org.mitk.views.modules");

  auto advisor = new QmitkM2aiaAppWorkbenchWindowAdvisor(this, configurer);
  advisor->ShowViewMenuItem(true);
  advisor->ShowNewWindowMenuItem(true);
  advisor->ShowClosePerspectiveMenuItem(false);
  advisor->SetPerspectiveExcludeList(perspExcludeList);
  advisor->SetViewExcludeList(viewExcludeList);
  advisor->ShowViewToolbar(false);  
  advisor->ShowPerspectiveToolbar(false);
  advisor->ShowVersionInfo(false);
  advisor->ShowMitkVersionInfo(false);
  advisor->ShowMemoryIndicator(true);
  advisor->SetProductName("M²aia Workbench");
  advisor->SetWindowIcon(":/M2aiaApplication/icon.png");

  return advisor;
}

QString QmitkM2aiaAppWorkbenchAdvisor::GetInitialWindowPerspectiveId()
{
  return "org.mitk.m2.application.perspectives.spectrumimaging";
}
