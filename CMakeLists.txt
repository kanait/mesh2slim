﻿if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "DO NOT BUILD in-tree.")
endif()

cmake_minimum_required ( VERSION 3.8 )
project( mesh2slim )

SET(CMAKE_C_COMPILER_WORKS 1 CACHE INTERNAL "")
SET(CMAKE_CXX_COMPILER_WORKS 1 CACHE INTERNAL "")

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

target_include_directories( mesh2slim
                            PUBLIC
                            ${PROJECT_SOURCE_DIR}/external/vecmath-cpp
                            ${PROJECT_SOURCE_DIR}/external/eigen
)

target_compile_options( mesh2slim PUBLIC -O3 -Wall)
target_compile_features( mesh2slim PUBLIC cxx_std_11 )

target_link_libraries( mesh2slim
                       LINK_PUBLIC
)

# Set StartUp Project (Option)
set_property( DIRECTORY PROPERTY VS_STARTUP_PROJECT "mesh2slim" )

