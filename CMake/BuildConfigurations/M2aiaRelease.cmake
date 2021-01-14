set(MITK_CONFIG_PACKAGES
  MatchPoint
  Qt5
  BLUEBERRY
)

set(MITK_CONFIG_PLUGINS
  org.mitk.gui.qt.datamanager
  org.mitk.gui.qt.stdmultiwidgeteditor
  org.mitk.planarfigure  
  org.mitk.gui.qt.imagecropper
  org.mitk.gui.qt.measurementtoolbox
  org.mitk.gui.qt.segmentation
  org.mitk.gui.qt.volumevisualization
  org.mitk.gui.qt.pointsetinteraction
  org.mitk.gui.qt.imagenavigator
  org.mitk.gui.qt.moviemaker
  org.mitk.gui.qt.properties
  org.mitk.gui.qt.viewnavigator
  org.mitk.gui.qt.renderwindowmanager
  org.mitk.gui.qt.multilabelsegmentation
  org.mitk.gui.qt.m2.application
  org.mitk.gui.qt.m2.Data
  org.mitk.gui.qt.m2.Ions
  org.mitk.gui.qt.m2.CombineImages
  org.mitk.gui.qt.m2.Spectrum
  org.mitk.gui.qt.m2.PeakPicking
  org.mitk.gui.qt.m2.ImzMLExport
  org.mitk.gui.qt.m2.Reconstruction3D
  org.mitk.gui.qt.m2.Position
  org.mitk.gui.qt.m2.Manipulator
  org.mitk.gui.qt.m2.OpticalImageRegistration
  org.mitk.m2.core.helper
)

set(MITK_VTK_DEBUG_LEAKS OFF CACHE BOOL "Enable VTK Debug Leaks" FORCE)

find_package(Doxygen REQUIRED)

# Ensure that the in-application help can be build
set(BLUEBERRY_QT_HELP_REQUIRED ON CACHE BOOL "Required Qt help documentation in plug-ins" FORCE)
