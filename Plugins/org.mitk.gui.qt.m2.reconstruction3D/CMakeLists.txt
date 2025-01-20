project(org_mitk_gui_qt_m2_reconstruction3D)

mitk_create_plugin(
  EXPORT_DIRECTIVE m2Reconstruction3D_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt M2aiaCore MitkElastix
)
