##Script to find CombineHarvester (mainly copied from https://gitlab.cern.ch/cms-desy-top/TopAnalysis/-/blob/master/Configuration/analysis/cmake/FindCommonClassifier.cmake)
SET(CMSSW_DIR "$ENV{CMSSW_BASE}")

FILE(GLOB_RECURSE localLibraries FOLLOW_SYMLINKS ${CMSSW_DIR}/lib/.poisonededmplugincache )
LIST(GET localLibraries 0 anyLibrary )
GET_FILENAME_COMPONENT(CMSSW_LIBRARIES_PATH ${anyLibrary} PATH)


FIND_LIBRARY( CombineTools_LIBRARY
  NAMES CombineHarvesterCombineTools
  PATHS ${CMSSW_LIBRARIES_PATH}
  NO_DEFAULT_PATH )
MARK_AS_ADVANCED(CombineTools_LIBRARY)

SET(CombineTools_INCLUDE_DIR ${CMSSW_DIR}/src)
MARK_AS_ADVANCED(CombineTools_INCLUDE_DIR)

IF ( CombineTools_INCLUDE_DIR AND CombineTools_LIBRARY )
  SET(CombineTools_FOUND TRUE)
  SET(CombineTools_LIBRARIES ${CombineTools_LIBRARY} )
  SET(CombineTools_INCLUDE_DIRS ${CombineTools_INCLUDE_DIR} )
ENDIF ( CombineTools_INCLUDE_DIR AND CombineTools_LIBRARY )


MESSAGE(STATUS "--------------- CombineTools SUMMARY ----------------")
MESSAGE(STATUS "CMSSW_DIR:              " ${CMSSW_DIR} )
MESSAGE(STATUS "CMSSW_LIBRARIES_PATH:   " ${CMSSW_LIBRARIES_PATH} )
MESSAGE(STATUS "CombineTools_INCLUDE_DIRS: " ${CombineTools_INCLUDE_DIRS} )
MESSAGE(STATUS "CombineTools_LIBRARIES:    " "${CombineTools_LIBRARIES}" )
MESSAGE(STATUS "CombineTools_DEFINITIONS:  " ${CombineTools_DEFINITIONS} )
MESSAGE(STATUS "---------------------------------------------")

IF ( CombineTools_FIND_REQUIRED AND NOT CombineTools_FOUND )
  MESSAGE( FATAL_ERROR "Cannot find CombineTools!!" )
ENDIF ( CombineTools_FIND_REQUIRED AND NOT CombineTools_FOUND )

link_directories("${CombineTools_LIBRARIES}")
include_directories("${CombineTools_INCLUDE_DIRS}")
