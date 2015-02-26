# CMakeFiles.txt for volregrid
#
# Andrew Janke - a.janke@gmail.com



CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(VOLREGRID)

#INCLUDE( ${LIBMINC_USE_FILE} )
ADD_DEFINITIONS(-DHAVE_CONFIG_H)

SET(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")


# check for prereqs
INCLUDE(CheckFunctionExists)
#CHECK_FUNCTION_EXISTS(mkstemp  HAVE_MKSTEMP)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_SOURCE_DIR}/config.h @ONLY)  

FIND_PACKAGE(FLEX)
#FIND_PACKAGE(LIBMINC REQUIRED)

FLEX_TARGET(lex lex.l lex.c )


ADD_EXECUTABLE(volregrid
                volregrid.c
                arb_path_io.h
                arb_path_io.c
                minc_support.h
                minc_support.c
                ${FLEX_lex_OUTPUTS}
               )

TARGET_INCLUDE_DIRECTORIES(volregrid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} /opt/minc/include )
TARGET_LINK_LIBRARIES(volregrid ${FLEX_LIBRARIES} m /opt/minc/lib )
                                                                                                                                                                    
INSTALL( TARGETS volregrid  DESTINATION bin)
