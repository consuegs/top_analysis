set(MODULES
  distributions
  # ~unfolding
  # ~unfolding2D
  # ~plot2D_diffMET
  # ~plot_diffMET
  # ~check_DM
  plot_distributions
  plot_distributions_combined
  # ~templateFit
  # ~compare_ttMC
  # ~test_MT2
  # ~plot2D_sensitivity
  # ~plot2D_sensitivity_MT
  # ~binningUnfolding
  # ~binningUnfolding_MT
  # ~scaleFactorUnfolding
  # ~plot_scaleFactorUnfolding
  # ~plot_response
  # ~compare_pseudoGen
  # ~diff_dPhiResolution
  # ~plot_dPhiResolution
  TUnfold_binning
  TUnfold_unfolding
  TUnfold_plotting
  TUnfold_combineYears
  # ~TUnfold_minimalTest
  # ~RooUnfold_test
  # ~plot_metSF
  # ~TUnfold_binning_fakeBin
  # ~TUnfold_unfolding_fakeBin
  # ~TUnfold_plotting_fakeBin
  # ~diff_fakes
  # ~plot_diff_fakes
  # ~diff_MET
  plot_metStudies
  # ~reweightingSF
  # ~test_TMVA
  # ~distributions_Philipp
  # ~plot_distributions_Philipp
  # ~Routin
  # ~test_BJetRegression
  triggerEff
  plot_triggerEff
  # ~dnn_Classification
  # ~dnn_Classification_fourVec
  # ~binningUnfolding_DNN
  # ~binningUnfolding_DNN_multi
  # ~binningUnfolding_DNN_cppflow
  # ~binningUnfolding_DNN_cppflow_2D
  # ~plot_DNNresolution
  # ~applyDNN
  # ~check_leptonVeto
  # ~test_tensorflow
  # ~plot_leptonVeto
  # ~plot_yearComparisons
  compare_powheg_amcatnlo
  binningUnfolding1D_Fabian_OA
  binningUnfolding1D_Fabian_OA_dPhi
  binningUnfolding2D_Fabian_OA
  bTagEff
  plot_bTagEff
  # ~distributions_KITsync
  # ~plot_yearComparisonsSebastian
  cr_envelope
  analyze_minTrees
  # ~testHarvester
  # ~testHarvester_ttbar
  # ~testHarvester_ttbar_chan
  # ~testHarvester_ttbar_2D
  # ~testONNX
  datacard_harvester
  datacard_harvester_2D
  )

foreach(MODULE ${MODULES})
  add_library(${MODULE} SHARED ${MODULE}.cpp)
endforeach()
