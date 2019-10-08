set(MODULES
  distributions
  unfolding
  plot2D_diffMET
  plot_diffMET
  check_DM
  plot_distributions
  templateFit
  compare_ttMC
  test_MT2
  plot2D_sensitivity
  binningUnfolding
  )

foreach(MODULE ${MODULES})
  add_library(${MODULE} SHARED ${MODULE}.cpp)
endforeach()
