if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "DO NOT BUILD in-tree.")
endif()

cmake_minimum_required ( VERSION 3.8 )
project( mesh2slim )

SET(CMAKE_C_COMPILER_WORKS 1 CACHE INTERNAL "")
SET(CMAKE_CXX_COMPILER_WORKS 1 CACHE INTERNAL "")

option(USE_SYSTEM_EIGEN "Prefer Eigen from system (Homebrew, etc.)" ON)

add_executable( mesh2slim
                BasisFunction.hxx
                MeshR.hxx
                MeshRApp.hxx
                MeshRDG.hxx
                PIEM.cxx
                PIEM.hxx
                PIEMHandler.hxx
                RIO.hxx
                SMFRIO.hxx
                SlimBall.hxx
                SlimSimplify.cxx
                SlimSimplify.hxx
                SlimTree.hxx
                SlimTreeApp.hxx
                SlimTreeIO.hxx
                envDep.h
                mesh2slim.cxx
                mydef.h
                pq.h
                strutil.h
                timer.hxx
                tokenizer.h
)

if (USE_SYSTEM_EIGEN)
  find_package(Eigen3 3.3 QUIET NO_MODULE)
  if (NOT TARGET Eigen3::Eigen)
    # Homebrew's eigen formula may not ship an Eigen3Config.cmake.
    # Fall back to locating headers (e.g. /opt/homebrew/include/eigen3/Eigen/Core).
    find_path(EIGEN3_INCLUDE_DIR
              NAMES Eigen/Core
              PATH_SUFFIXES eigen3 include/eigen3)
  endif()
endif()

target_include_directories( mesh2slim
                            PUBLIC
                            ${PROJECT_SOURCE_DIR}/external/vecmath-cpp
)

if (TARGET Eigen3::Eigen)
  target_link_libraries(mesh2slim PUBLIC Eigen3::Eigen)
elseif (EIGEN3_INCLUDE_DIR)
  message(STATUS "Using system Eigen headers at: ${EIGEN3_INCLUDE_DIR}")
  target_include_directories(mesh2slim PUBLIC ${EIGEN3_INCLUDE_DIR})
else()
  message(STATUS "Eigen3 not found via find_package(); falling back to external/eigen")
  target_include_directories(mesh2slim PUBLIC ${PROJECT_SOURCE_DIR}/external/eigen)
endif()

target_compile_options( mesh2slim PUBLIC -O3 -Wall)
target_compile_features( mesh2slim PUBLIC cxx_std_11 )

target_link_libraries( mesh2slim
                       LINK_PUBLIC
)

# Set StartUp Project (Option)
set_property( DIRECTORY PROPERTY VS_STARTUP_PROJECT "mesh2slim" )

