########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( mymarlin )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( MarlinUtil 1.4 REQUIRED )
FIND_PACKAGE( MarlinKinfit 0.0 REQUIRED )
FIND_PACKAGE( MarlinTrk )
FIND_PACKAGE( GSL REQUIRED )
FIND_PACKAGE( ROOT REQUIRED )
FIND_PACKAGE( DD4hep COMPONENTS DDRec )
FIND_PACKAGE( aidaTT )
FIND_PACKAGE( GBL REQUIRED)

OPTION( MARLINRECO_AIDA "Set to ON to build MarlinReco with AIDA" ON )
IF( MARLINRECO_AIDA )
    FIND_PACKAGE( AIDA )
    IF( NOT AIDA_FOUND )
        MESSAGE( SEND_ERROR "AIDA not found. Please set MARLINRECO_AIDA to OFF or AIDA_DIR=/path/to/AIDA" )
    ENDIF()
ENDIF()

SET( ROOT_DEFINITIONS "-DUSEROOT -DUSE_ROOT -DMARLIN_USE_ROOT" )

FOREACH( pkg Marlin MarlinUtil MarlinKinfit MarlinTrk GSL AIDA ROOT aidaTT GBL)
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( SYSTEM ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()




### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################

# definitions to pass to the compiler
# ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
# ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ./include )
#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

