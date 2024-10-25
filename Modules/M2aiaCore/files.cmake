set(H_FILES 
  include/m2CoreCommon.h
  include/m2Process.hpp 
  include/m2SpectrumImageHelper.h
  include/m2SpectrumImageStack.h
  include/m2SpectrumImageDataInteractor.h
  
  include/m2ImzMLParser.h
  include/m2Timer.h
  include/m2SpectrumInfo.h
  include/m2ImzMLImageIO.h
  include/m2ImzMLEngine.h
  include/m2TestFixture.h
  
  # include/m2ElxUtil.h
  # include/m2ElxRegistrationHelper.h
  include/m2SubdivideImage2DFilter.h
  include/m2ShiftMapImageFilter.h
  
  include/npy/npy.hpp
  
  include/m2SpectrumImage.h
  include/m2ISpectrumImageDataAccess.h
  include/m2ISpectrumImageSource.h
  
  include/m2ImzMLSpectrumImage.h
  include/m2ImzMLSpectrumImageSource.hpp
  
  include/m2FsmSpectrumImage.h

  include/m2IntervalVector.h

  include/m2DataNodePredicates.h

  include/signal/m2Baseline.h
  include/signal/m2EstimateFwhm.h
  include/signal/m2MedianAbsoluteDeviation.h
  include/signal/m2Morphology.h
  include/signal/m2Normalization.h
  include/signal/m2Binning.h
  include/signal/m2PeakDetection.h
  include/signal/m2Pooling.h
  include/signal/m2RunningMedian.h
  include/signal/m2SignalCommon.h
  include/signal/m2Smoothing.h
  include/signal/m2Transformer.h

)

set(CPP_FILES
  m2SpectrumImage.cpp
  m2SpectrumImageHelper.cpp
  m2SpectrumImageStack.cpp
  m2CoreObjectFactory.cpp 
  m2ImzMLSpectrumImage.cpp
  m2FsmSpectrumImage.cpp
  m2SubdivideImage2DFilter.cpp
  m2SpectrumImageDataInteractor.cpp
  m2IntervalVector.cpp
  m2ShiftMapImageFilter.cpp
  
  IO/m2ImzMLParser.cpp
  IO/m2ImzMLImageIO.cpp
  IO/m2ImzMLEngine.cpp
  IO/m2PythonWrapper.cpp
)

set(RESOURCE_FILES
  Interactions/SpectrumDataInteractor.xml
  Interactions/SpectrumDataInteractorConfig.xml
)
