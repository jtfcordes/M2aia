project(org_mitk_gui_qt_m2_application)

mitk_create_plugin(
  SUBPROJECTS MITK-M2aia
  EXPORT_DIRECTIVE M2AIA_APP_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgets MitkSceneSerialization
  PACKAGE_DEPENDS Qt6|WebEngineWidgets
)
