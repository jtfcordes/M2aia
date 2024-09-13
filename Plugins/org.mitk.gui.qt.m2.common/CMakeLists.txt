project(org_mitk_gui_qt_m2_common)

mitk_create_plugin(
  EXPORT_DIRECTIVE MITK_M2_CORE_HELPER_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgets MitkQtWidgetsExt M2aiaCore
  PACKAGE_DEPENDS Qt6|UiTools
)
