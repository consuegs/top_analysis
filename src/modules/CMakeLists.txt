set(MODULES
  distributions
  unfolding
  unfolding2D
  plot2D_diffMET
  plot_diffMET
  check_DM
  plot_distributions
  templateFit
  compare_ttMC
  test_MT2
  plot2D_sensitivity
  plot2D_sensitivity_MT
  binningUnfolding
  binningUnfolding_MT
  scaleFactorUnfolding
  plot_scaleFactorUnfolding
  plot_response
  compare_pseudoGen
  diff_dPhiResolution
  plot_dPhiResolution
  TUnfold_binning
  TUnfold_unfolding
  TUnfold_plotting
  TUnfold_minimalTest
  RooUnfold_test
  plot_metSF
  )

foreach(MODULE ${MODULES})
  add_library(${MODULE} SHARED ${MODULE}.cpp)
endforeach()
