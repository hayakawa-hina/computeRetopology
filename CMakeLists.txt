PROJECT(patchgen_demo)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE_DIRECTORIES(./src)
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/lib)

# Windows compiler settings
IF(WIN32)
  ADD_DEFINITIONS(
    -D__STDC_LIMIT_MACROS
    -DEIGEN_DONT_ALIGN_STATICALLY
    -DGLEW_STATIC
    -DFREEGLUT_STATIC
  )
  SET(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   /MTd")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT ")
  SET(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   /MTd")
  SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT ")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
ENDIF()

# Mac/Linux compiler settings
IF(UNIX)
  INCLUDE_DIRECTORIES(
    /usr/local/include
    /opt/local/include
  )
  LINK_DIRECTORIES(
    /usr/local/lib
    /opt/local/lib
  )
  IF(APPLE)
    SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 -stdlib=libc++ -ObjC++ -framework OpenGL -framework GLUT -framework Cocoa")
    ADD_DEFINITIONS(-D_MACOSX)
  ELSE()
    SET(CMAKE_CXX_FLAGS "-std=c++11 -O3")
    ADD_DEFINITIONS(-D_UNIX)
  ENDIF()
ENDIF()

# AntTweakBar
ADD_DEFINITIONS(
  -DTW_NO_LIB_PRAGMA
  -DTW_STATIC
)
INCLUDE_DIRECTORIES(./src/AntTweakBar)
FILE(GLOB_RECURSE anttweakbar_src
  ./src/AntTweakBar/*.c
  ./src/AntTweakBar/*.cpp
)
ADD_LIBRARY(AntTweakBar ${anttweakbar_src})

# lp_solve_5.5
ADD_DEFINITIONS(
  -DLoadInverseLib=FALSE
  -DLoadLanguageLib=FALSE
  -DLoadableBlasLib=0
  -DINVERSE_ACTIVE=INVERSE_LUSOL
  -DRoleIsExternalInvEngine
  -DNOFORTIFY
  -DCHECK_SOLUTION
  -DYY_NEVER_INTERACTIVE
  -DPARSER_LP
  -DLP_MAXLINELEN=0
)
INCLUDE_DIRECTORIES(
  ./src/lp_solve_5.5
  ./src/lp_solve_5.5/bfp
  ./src/lp_solve_5.5/bfp/bfp_LUSOL
  ./src/lp_solve_5.5/bfp/bfp_LUSOL/LUSOL
  ./src/lp_solve_5.5/colamd
  ./src/lp_solve_5.5/shared
)
FILE(GLOB lp_solve_src ./src/lp_solve_5.5/*.h ./src/lp_solve_5.5/*.c)
SET(lp_solve_bfp_src
  ./src/lp_solve_5.5/bfp/bfp_LUSOL/lp_LUSOL.c
  ./src/lp_solve_5.5/bfp/bfp_LUSOL/LUSOL/lusol.c
)
SET(lp_solve_colamd_src
  ./src/lp_solve_5.5/colamd/colamd.h
  ./src/lp_solve_5.5/colamd/colamd.c
)
FILE(GLOB lp_solve_shared_src ./src/lp_solve_5.5/shared/*.h ./src/lp_solve_5.5/shared/*.c)
SOURCE_GROUP("lp_solve_5.5" FILES ${lp_solve_src})
SOURCE_GROUP("lp_solve_5.5\\bfp" FILES ${lp_solve_bfp_src})
SOURCE_GROUP("lp_solve_5.5\\colamd" FILES ${lp_solve_colamd_src})
SOURCE_GROUP("lp_solve_5.5\\shared" FILES ${lp_solve_shared_src})
ADD_LIBRARY(lp_solve
  ${lp_solve_src}
  ${lp_solve_bfp_src}
  ${lp_solve_colamd_src}
  ${lp_solve_shared_src}
)

# OpenMesh
ADD_DEFINITIONS(-DOM_STATIC_BUILD)
IF(WIN32)
  ADD_DEFINITIONS(
    -D_USE_MATH_DEFINES
    -DNOMINMAX
  )
ENDIF(WIN32)
FILE(GLOB_RECURSE openmesh_src ./src/OpenMesh/*.cc)
ADD_LIBRARY(OpenMesh ${openmesh_src})

# patchgen
FILE(GLOB patchgen_src ./src/patchgen/*.h ./src/patchgen/*.cpp)
ADD_LIBRARY(patchgen ${patchgen_src})

# patchgen_demo
FILE(GLOB patchgen_demo_src ./src/patchgen_demo/*.h ./src/patchgen_demo/*.cpp)
ADD_EXECUTABLE(patchgen_demo
  ${patchgen_demo_src})
IF(WIN32)
  TARGET_LINK_LIBRARIES(patchgen_demo
    patchgen
    AntTweakBar
    lp_solve
    OpenMesh
    OpenGL32
    GLU32
    freeglut_static
    glew32
  )
ELSEIF(APPLE)
    FIND_LIBRARY(Cocoa_LIB Cocoa)
    FIND_LIBRARY(GLUT_LIB GLUT)
    FIND_LIBRARY(GLEW_LIB GLEW)
    FIND_LIBRARY(OpenGL_LIB OpenGL)
    MARK_AS_ADVANCED(Cocoa_LIB GLUT_LIB GLEW_LIB OpenGL_LIB)
  TARGET_LINK_LIBRARIES(patchgen_demo
    patchgen
    AntTweakBar
    lp_solve
    OpenMesh
    patchgen
    ${Cocoa_LIB}
    ${GLUT_LIB}
    ${GLEW_LIB}
    ${OpenGL_LIB}
  )
ELSEIF(UNIX)
  TARGET_LINK_LIBRARIES(patchgen_demo
    patchgen
    AntTweakBar
    lp_solve
    OpenMesh
    glut
    GLEW
    GLU
    GL
    X11
    pthread
  )
ENDIF()
