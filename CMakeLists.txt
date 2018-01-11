cmake_minimum_required(VERSION 2.8.9)
PROJECT(eltopo)

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
ADD_DEFINITIONS(-D__LITTLE_ENDIAN__ -DUSE_FORTRAN_BLAS)

INCLUDE_DIRECTORIES(eltopo3d common common/newsparse common/meshes common/tunicate)

FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)



ADD_LIBRARY(eltopo_common STATIC ${ELTOPO_COMMON_SRC})
ADD_LIBRARY(eltopo STATIC ${ELTOPO_SRC})
TARGET_LINK_LIBRARIES(eltopo eltopo_common ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})






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
TARGET_LINK_LIBRARIES(talpa eltopo ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})

