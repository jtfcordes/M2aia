project(org_mitk_gui_qt_m2_spectrum)

mitk_create_plugin(
  EXPORT_DIRECTIVE M2_SPECTRUM_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt M2aiaCore
  PACKAGE_DEPENDS Qt6|Charts Boost
)
