cmake_minimum_required(VERSION 2.8.9)
PROJECT(eltopo)



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -Wextra -Wno-reorder -Wno-unknown-pragmas -march=native" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")


OPTION(GUI_ENABLED "USE GUI" OFF)
OPTION(SOLVERS_ENABLED "USE SOLVERS" ON)

SET(SOLVERS_SRC
    common/newsparse/sparse_matrix.cpp
    common/newsparse/krylov_solvers.cpp
    )

# Common
SET(ELTOPO_COMMON_SRC
    common/bfstream.cpp
    common/collisionqueries.cpp
    common/cubic_ccd_wrapper.cpp
    common/fileio.cpp
    common/levelset.cpp
    common/makelevelset2.cpp
    common/makelevelset3.cpp
    common/marching_tiles_hires.cpp
    common/marching_triangles.cpp
    common/root_parity_ccd_wrapper.cpp
    common/runstats.cpp
    common/wallclocktime.cpp
    common/tunicate/expansion.cpp
    common/tunicate/intersection.cpp
    common/tunicate/interval.cpp
    common/tunicate/neg.cpp
    common/tunicate/orientation.cpp
    common/tunicate/rootparitycollisiontest.cpp
    common/newsparse/dense_matrix.cpp
    common/newsparse/krylov_solvers.cpp
    common/newsparse/sparse_matrix.cpp
    )
IF(SOLVERS_ENABLED)
    FIND_PACKAGE(BLAS REQUIRED)
    FIND_PACKAGE(LAPACK REQUIRED)

    include_directories(${BLAS_INCLUDE_DIR})
    include_directories(${LAPACK_INCLUDE_DIR})
    SET(ELTOPO_COMMON_SRC ${ELTOPO_COMMON_SRC} ${SOLVERS_SRC})
ENDIF(SOLVERS_ENABLED)
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
    eltopo3d/nondestructivetrimesh.cpp 
    eltopo3d/subdivisionscheme.cpp 
    eltopo3d/surftrack.cpp
    eltopo3d/trianglequality.cpp
    )



IF(GUI_ENABLED)
    SET(OpenGL_GL_PREFERENCE GLVND)
    FIND_PACKAGE(OpenGL REQUIRED)
    FIND_PACKAGE(GLUT REQUIRED)
    INCLUDE_DIRECTORIES( ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} )
    ADD_DEFINITIONS(-DUSE_GUI)

    SET(ELTOPO_COMMON_SRC ${ELTOPO_COMMON_SRC} common/gluvi.cpp)
    SET(ELTOPO_SRC 
        ${ELTOPO_SRC}
        eltopo3d/meshrenderer.cpp 
        )

ELSE(GUI_ENABLED)
    ADD_DEFINITIONS(-DNO_GUI)
    SET(OPENGL_LIBRARIES )
    SET(GLUT_LIBRARY )
ENDIF(GUI_ENABLED)


find_package(Threads REQUIRED)


ADD_LIBRARY(eltopo_common STATIC ${ELTOPO_COMMON_SRC})
TARGET_INCLUDE_DIRECTORIES(eltopo_common PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/common 
    ${CMAKE_CURRENT_SOURCE_DIR}/common/newsparse 
    ${CMAKE_CURRENT_SOURCE_DIR}/common/meshes 
    ${CMAKE_CURRENT_SOURCE_DIR}/common/tunicate)


if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    target_compile_definitions(eltopo_common PUBLIC USE_FORTRAN_BLAS)
    target_compile_definitions(eltopo_common PUBLIC __LITTLE_ENDIAN__)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flax-vector-conversions")
endif()


ADD_LIBRARY(eltopo STATIC ${ELTOPO_SRC})

TARGET_INCLUDE_DIRECTORIES(eltopo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/eltopo3d)
TARGET_INCLUDE_DIRECTORIES(eltopo INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
TARGET_LINK_LIBRARIES(eltopo eltopo_common ${GLUT_LIBRARY} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(GUI_ENABLED)
    TARGET_LINK_LIBRARIES(eltopo OpenGL::GL OpenGL::GLU)
endif()


if(CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops")
endif()

# We need C++11. Put this directive after CGAL's include.
if (APPLE)
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()
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

ADD_EXECUTABLE(talpa EXCLUDE_FROM_ALL ${TALPA_SRC})
TARGET_INCLUDE_DIRECTORIES(talpa PUBLIC talpa talpa/drivers talpa/curlnoise)
TARGET_LINK_LIBRARIES(talpa eltopo Threads::Threads)

