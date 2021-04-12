set(HEADERS
   cppflow/context.h
   cppflow/cppflow.h
   cppflow/datatype.h
   cppflow/defer.h
   cppflow/model.h
   cppflow/ops.h
   cppflow/raw_ops.h
   cppflow/tensor.h
  )


root_generate_dictionary(cppflow_dict ${HEADERS} LINKDEF cppflow/cppflow_LinkDef.h)

add_library(cppflow SHARED cppflow.cpp cppflow_dict)
target_link_libraries (cppflow "${TENSORFLOW_LIB}")
target_link_libraries (cppflow "${Protobuf_LIB}")
target_link_libraries (cppflow "${Giflib_LIB}")
add_custom_command(TARGET cppflow POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/src/cppflow/cppflow_dict_rdict.pcm ${LIBRARY_OUTPUT_PATH}
  )

