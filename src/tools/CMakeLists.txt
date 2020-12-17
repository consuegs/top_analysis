set(SRCS
  io.cpp
  util.cpp
  physics.cpp
  gfx.cpp
  hist.cpp
  efficiency.cpp
  MT2Functor.cc
  selection.cpp
  )

add_library(tools SHARED ${SRCS})

