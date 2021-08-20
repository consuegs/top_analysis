set(SRCS
   JetCorrectionUncertainty.cc
   JetCorrectorParameters.cc
   SimpleJetCorrectionUncertainty.cc
   Utilities.cc
   JetResolution.cc
   JetResolutionObject.cc
  )

add_library(ext SHARED ${SRCS})

