cmake_minimum_required(VERSION 2.8.9)
PROJECT(eltopo)


if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -g -Wall -Wextra -Wno-reorder -Wno-unknown-pragmas -march=native" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")


OPTION(GUI_ENABLED "USE GUI" OFF)

# Common
SET(ELTOPO_COMMON_SRC common/wallclocktime.cpp
    common/newsparse/sparse_matrix.cpp
    common/newsparse/krylov_solvers.cpp

    common/tunicate/expansion.cpp common/tunicate/intersection.cpp common/tunicate/neg.cpp
    common/tunicate/orientation.cpp

    common/root_parity_ccd_wrapper.cpp common/cubic_ccd_wrapper.cpp common/collisionqueries.cpp
    common/bfstream.cpp common/makelevelset3.cpp common/marching_tiles_hires.cpp
    common/runstats.cpp 
    )
# Source files
SET(ELTOPO_SRC 
    eltopo3d/accelerationgrid.cpp 
    eltopo3d/broadphasegrid.cpp 
    eltopo3d/collisionpipeline.cpp
    eltopo3d/dynamicsurface.cpp 
    eltopo3d/edgecollapser.cpp 
    eltopo3d/edgeflipper.cpp 
    eltopo3d/edgesplitter.cpp
    eltopo3d/eltopo.cpp 
    eltopo3d/impactzonesolver.cpp 
    eltopo3d/meshmerger.cpp 
    eltopo3d/meshpincher.cpp 
    eltopo3d/meshsmoother.cpp
    eltopo3d/meshrenderer.cpp 
    eltopo3d/nondestructivetrimesh.cpp 
    eltopo3d/subdivisionscheme.cpp 
    eltopo3d/surftrack.cpp
    eltopo3d/trianglequality.cpp
    )



IF(GUI_ENABLED)
FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(GLUT REQUIRED)
INCLUDE_DIRECTORIES( ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} )
ADD_DEFINITIONS(-DUSE_GUI)
ELSE(GUI_ENABLED)
ADD_DEFINITIONS(-DNO_GUI)
SET(OPENGL_LIBRARIES )
SET(GLUT_LIBRARY )
ENDIF(GUI_ENABLED)

INCLUDE_DIRECTORIES(eltopo3d 
    ${CMAKE_CURRENT_SOURCE_DIR}/common 
    ${CMAKE_CURRENT_SOURCE_DIR}/common/newsparse 
    ${CMAKE_CURRENT_SOURCE_DIR}/common/meshes 
    ${CMAKE_CURRENT_SOURCE_DIR}/common/tunicate)

FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

include_directories(${BLAS_INCLUDE_DIR})
include_directories(${LAPACK_INCLUDE_DIR})

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_definitions(-DUSE_FORTRAN_BLAS)
  add_definitions(-D__LITTLE_ENDIAN__)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flax-vector-conversions")
endif()




set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops")

ADD_LIBRARY(eltopo_common STATIC ${ELTOPO_COMMON_SRC})
ADD_LIBRARY(eltopo STATIC ${ELTOPO_SRC})
TARGET_LINK_LIBRARIES(eltopo eltopo_common ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# We need C++11. Put this directive after CGAL's include.
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()



SET(TALPA_SRC 
    talpa/iomesh.cpp 
    talpa/lexer.cpp 
    talpa/newparser.cpp 
    talpa/geometryinit.cpp 
    talpa/scriptinit.cpp 
    # drivers
    talpa/curlnoise/noise.cpp 
    talpa/drivers/meancurvature.cpp 
    talpa/drivers/faceoff.cpp 
    talpa/drivers/normaldriver.cpp 
    talpa/drivers/sisccurlnoisedriver.cpp
    talpa/main.cpp
    )

ADD_EXECUTABLE(talpa ${TALPA_SRC})
TARGET_INCLUDE_DIRECTORIES(talpa PUBLIC talpa talpa/drivers talpa/curlnoise)
TARGET_LINK_LIBRARIES(talpa eltopo)



SET(PYBIND11_CPP_STANDARD -std=c++1z)

ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/extern/pybind11)


#ADD_LIBRARY(pyeltopo MODULE pyeltopo/main.cpp
#    pyeltopo/eltopo.cpp)
PYBIND11_ADD_MODULE(pyeltopo pyeltopo/main.cpp
    pyeltopo/eltopo.cpp)

TARGET_INCLUDE_DIRECTORIES(pyeltopo PRIVATE .)
TARGET_INCLUDE_DIRECTORIES(pyeltopo PRIVATE /usr/include/eigen3)
TARGET_INCLUDE_DIRECTORIES(pyeltopo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/extern/core/include)

target_link_libraries(pyeltopo PRIVATE eltopo)
#target_link_libraries(pyeltopo PRIVATE pybind11::module eltopo)
#set_target_properties(pyeltopo PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
#                                         SUFFIX "${PYTHON_MODULE_EXTENSION}")

MESSAGE(STATUS "Pyexec: ${PYTHON_EXECUTABLE}")

EXECUTE_PROCESS(
    COMMAND "${PYTHON_EXECUTABLE}" -m site --user-site
    OUTPUT_VARIABLE PYTHON_SITE
    OUTPUT_STRIP_TRAILING_WHITESPACE)
MESSAGE(STATUS "Pysites: ${PYTHON_SITE}")

MESSAGE(STATUS "TARGETS pyeltopo LIBRARY DESTINATION ${PYTHON_SITE}")
INSTALL(TARGETS pyeltopo LIBRARY DESTINATION ${PYTHON_SITE})
