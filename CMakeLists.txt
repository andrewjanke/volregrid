# CMakeFiles.txt for volregrid
#
# Andrew Janke - a.janke@gmail.com
#
# do this:
# export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/opt/minc/lib
#
# to ensure that CMake finds /opt/minc/lib/LIBMINCConfig.cmake

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


PROJECT(volregrid)

# get current version from git tag
EXECUTE_PROCESS(COMMAND git describe
   OUTPUT_VARIABLE VERSION_STRING
   OUTPUT_STRIP_TRAILING_WHITESPACE)

ADD_DEFINITIONS(-DPACKAGE_NAME="${CMAKE_PROJECT_NAME}")
ADD_DEFINITIONS(-DPACKAGE_VERSION="${VERSION_STRING}")
ADD_DEFINITIONS(-DPACKAGE_BUGREPORT="a.janke@gmail.com")

# load any modules used in this project
SET(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")

# find the pre-requisites
FIND_PACKAGE(LIBMINC REQUIRED)
FIND_PACKAGE(FLEX REQUIRED)
FIND_PACKAGE(GSL REQUIRED)

# set compile options
INCLUDE( ${LIBMINC_USE_FILE} )
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
LINK_LIBRARIES(${LIBMINC_LIBRARIES} ${GSL_LIBRARIES})

# build targets
FLEX_TARGET(lex lex.l ${CMAKE_CURRENT_BINARY_DIR}/lex.c )

ADD_EXECUTABLE(volregrid
                volregrid.c
                arb_path_io.h
                arb_path_io.c
                minc_support.h
                minc_support.c
                ${FLEX_lex_OUTPUTS}
               )

# what and where to install
INSTALL( TARGETS volregrid  DESTINATION bin)


