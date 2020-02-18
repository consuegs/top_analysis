set(HEADERS
   RooUnfold/RooUnfold.h
   RooUnfold/RooUnfoldBayes.h
   RooUnfold/RooUnfoldBinByBin.h
   RooUnfold/RooUnfoldErrors.h
   RooUnfold/RooUnfoldIds.h
   RooUnfold/RooUnfoldInvert.h
   RooUnfold/RooUnfoldParms.h
   RooUnfold/RooUnfoldResponse.h
   RooUnfold/RooUnfoldSvd.h
   RooUnfold/RooUnfoldTUnfold.h
   RooUnfold/TSVDUnfold_local.h
  )

root_generate_dictionary(RooUnfold_dict ${HEADERS} LINKDEF RooUnfold/RooUnfold_LinkDef.h)

set(SRCS
   RooUnfold.cxx
   RooUnfoldBayes.cxx
   RooUnfoldBinByBin.cxx
   RooUnfoldErrors.cxx
   RooUnfoldIds.cxx
   RooUnfoldInvert.cxx
   RooUnfoldParms.cxx
   RooUnfoldResponse.cxx
   RooUnfoldSvd.cxx
   RooUnfoldTUnfold.cxx
   TSVDUnfold.cxx
  )

add_library(RooUnfold SHARED ${SRCS} RooUnfold_dict)

add_custom_command(TARGET RooUnfold POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/src/RooUnfold/RooUnfold_dict_rdict.pcm ${LIBRARY_OUTPUT_PATH}
  )

