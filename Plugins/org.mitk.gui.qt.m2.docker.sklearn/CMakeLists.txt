project(org_mitk_gui_qt_m2_docker_sklearn)

mitk_create_plugin(
  EXPORT_DIRECTIVE PEAKLEARNING_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS PRIVATE MitkQtWidgetsExt M2aiaCore MitkDocker
) 

# if(DOCKER_EXECUTABLE)
#   add_custom_target(docker_sklearn
#     COMMAND ${DOCKER_EXECUTABLE} build --quiet  -f ${CMAKE_CURRENT_LIST_DIR}/docker/Dockerfile -t ghcr.io/m2aia/sklearn ${CMAKE_CURRENT_LIST_DIR}/docker
#   )
#   add_dependencies(org_mitk_gui_qt_m2_docker_sklearn docker_sklearn)
# endif()