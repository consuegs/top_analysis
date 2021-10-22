set(SRCS
   JetCorrectionUncertainty.cc
   JetCorrectorParameters.cc
   SimpleJetCorrectionUncertainty.cc
   Utilities.cc
   JetResolution.cc
   JetResolutionObject.cc
   JetResolutionObject.cc
   BTagCalibrationStandalone.cpp
  )

add_library(ext SHARED ${SRCS})

