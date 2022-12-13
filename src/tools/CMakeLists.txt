set(SRCS
  io.cpp
  util.cpp
  physics.cpp
  gfx.cpp
  hist.cpp
  efficiency.cpp
  MT2Functor.cc
  selection.cpp
  dataTrigger2016.cpp
  dataTrigger2017.cpp
  dataTrigger2018.cpp
  jetCorrections.cpp
  systematics.cpp
  bTagWeights.cpp
  leptonSF.cpp
  leptonCorrections.cpp
  triggerSF.cpp
  mcWeights.cpp
  minTreeReader.cpp
  dnnRegression.cpp
  tunfoldPlottingHelper.cpp
  distributionsPlottingHelper.cpp
  jetPileupIDSF.cpp
  )

add_library(tools SHARED ${SRCS})

