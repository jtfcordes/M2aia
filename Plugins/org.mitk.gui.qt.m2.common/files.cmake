set(SRC_CPP_FILES
  m2UIUtils.cpp
  # Qm2AssociatedFilesDialog.cpp  
)

set(INTERNAL_CPP_FILES
  org_mitk_gui_qt_m2_common_Activator.cpp
  m2BrowserPreferencesPage.cpp
  QmitkDataNodeConvertToRGBImageAction.cpp
  QmitkDataNodeConvertPixelTypeAction.cpp
  QmitkDataNodeExportComponentAction.cpp
  QmitkDataNodeReimportImageAction.cpp
  QmitkDataNodePlotColorAction.cpp
)

set(UI_FILES
  src/internal/m2PreferencePage.ui
  src/Qm2NameDialogControls.ui
  # src/Qm2AssociatedFilesDialog.ui
)

set(H_FILES
)

set(MOC_H_FILES
  src/internal/org_mitk_gui_qt_m2_common_Activator.h
  src/internal/m2BrowserPreferencesPage.h
  src/m2UIUtils.h
  src/Qm2NameDialog.h
  # src/Qm2AssociatedFilesDialog.h
  src/internal/QmitkDataNodeConvertToRGBImageAction.h
  src/internal/QmitkDataNodeConvertPixelTypeAction.h
  src/internal/QmitkDataNodeExportComponentAction.h
  src/internal/QmitkDataNodePlotColorAction.h
  src/internal/QmitkDataNodeReimportImageAction.h
)

# list of resource files which can be used by the plug-in
# system without loading the plug-ins shared library,
# for example the icon used in the menu and tabs for the
# plug-in views in the workbench
set(CACHED_RESOURCE_FILES
  plugin.xml
  resources/MAPRegData.svg
)

# list of Qt .qrc files which contain additional resources
# specific to this plugin
set(QRC_FILES
  resources/QmitkM2aiaCore.qrc
)

set(CPP_FILES )

foreach(file ${SRC_CPP_FILES})
  set(CPP_FILES ${CPP_FILES} src/${file})
endforeach(file ${SRC_CPP_FILES})

foreach(file ${INTERNAL_CPP_FILES})
  set(CPP_FILES ${CPP_FILES} src/internal/${file})
endforeach(file ${INTERNAL_CPP_FILES})

