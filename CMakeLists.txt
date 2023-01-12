cmake_minimum_required(VERSION 2.8)

project(ma_analysis)

# set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

### where to find the libraries
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
### where to find the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}")

### add ROOT
# ~find_package(ROOT REQUIRED TreePlayer PyROOT Unfold XMLParser)
# ~find_package(ROOT REQUIRED TreePlayer PyROOT Unfold XMLParser TMVA)
find_package(ROOT REQUIRED TreePlayer PyROOT Unfold XMLParser TMVA PyMVA)
include(RootUseFile)

### add cppflow (including tensorflow and dependencies)
find_library(TENSORFLOW_LIB tensorflow HINT /cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/tensorflow/2.1.0-ghbfee/lib/)
find_library(Protobuf_LIB protobuf HINT /cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/protobuf/3.11.3/lib)
find_library(Giflib_LIB gif HINT /cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/giflib/5.2.0/lib)
include_directories("/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/tensorflow/2.1.0-ghbfee/include/")
include_directories("/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/eigen/d812f411c3f9-ghbfee/include/eigen3")
link_directories("${TENSORFLOW_LIB}")
link_directories("${Protobuf_LIB}")
link_directories("${Giflib_LIB}")
include_directories("/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/protobuf/3.11.3/include")

### add ONNXRuntime
include_directories("${CMAKE_SOURCE_DIR}/src/ext/onnxruntime-linux-x64-1.11.0/include/")
link_directories("${CMAKE_SOURCE_DIR}/src/ext/onnxruntime-linux-x64-1.11.0/lib")

# change version: cmake --help-module FindPythonLibs
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

#~ set(BOOST_INCLUDEDIR "/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-ikhhed/include/")
#~ set(Boost_LIBRARY_DIR "/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-ikhhed/lib/")
#~ set(BOOST_INCLUDEDIR "/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/boost/1.63.0/include/")
#~ set(Boost_LIBRARY_DIR "/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/boost/1.63.0/lib/")
set(BOOST_INCLUDEDIR "/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/boost/1.67.0/include/")
set(Boost_LIBRARY_DIR "/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/boost/1.67.0/lib/")
find_package(Boost 1.67 COMPONENTS program_options REQUIRED)
include_directories(${BOOST_INCLUDEDIR})
include_directories(${Boost_LIBRARY_DIR})

# add CombineHarvester
FIND_PACKAGE( CombineHarvester )

### c++17
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
#~ set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# ~set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

### define standalone
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSTANDALONE")

### ONNXRuntime
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lonnxruntime")

# source /cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/valgrind/3.15.0/etc/profile.d/init.sh
# valgrind --num-callers=30 --keep-debuginfo=yes --error-limit=no --undef-value-errors=no --tool=memcheck --leak-check=full --suppressions=$ROOTSYS/etc/valgrind-root.supp ./run.x plot_distributions_combined > valgrind.txt 2>&1

# enable address sanatizer
# ~set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
# ~set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

## std::filesystem
link_libraries(stdc++fs)

# pass some definitions for usage in code
add_definitions(
  -DLIBRARY_OUTPUT_PATH="${LIBRARY_OUTPUT_PATH}/"
  -DCMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}/"
  )

# add the actual source directory
add_subdirectory(src)

# run target to execute "make run"
add_custom_target(run
  DEPENDS run.x tests
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND time ${EXECUTABLE_OUTPUT_PATH}/run.x
  )


