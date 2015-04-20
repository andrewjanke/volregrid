# CMakeFiles.txt for volregrid
#
# Andrew Janke - a.janke@gmail.com
#
# do this:
# export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:/opt/minc/lib
#
# to ensure that CMake finds /opt/minc/lib/LIBMINCConfig.cmake


CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(VOLREGRID)


SET(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h @ONLY)  

FIND_PACKAGE(LIBMINC REQUIRED)
FIND_PACKAGE(FLEX REQUIRED)
FIND_PACKAGE(GSL REQUIRED)
INCLUDE( ${LIBMINC_USE_FILE} )
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

FLEX_TARGET(lex lex.l ${CMAKE_CURRENT_BINARY_DIR}/lex.c )

LINK_LIBRARIES(${LIBMINC_LIBRARIES} ${GSL_LIBRARIES})

ADD_EXECUTABLE(volregrid
                volregrid.c
                arb_path_io.h
                arb_path_io.c
                minc_support.h
                minc_support.c
                ${FLEX_lex_OUTPUTS}
               )


INSTALL( TARGETS volregrid  DESTINATION bin)


