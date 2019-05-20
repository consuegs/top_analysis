set(MODULES
  distributions
  unfolding
  plot2D_diffMET
  plot_diffMET
  check_DM
  plot_distributions
  templateFit
  )

foreach(MODULE ${MODULES})
  add_library(${MODULE} SHARED ${MODULE}.cpp)
endforeach()
