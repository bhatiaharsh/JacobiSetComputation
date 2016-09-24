# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory
# Written by Harsh Bhatia (bhatia4@llnl.gov).
# CODE-701045.
# All rights reserved.

# This file is part of JacobiSetComputation v1.0.

# For details, see https://github.com/bhatiaharsh/JacobiSetComputation.
# For more details on the Licence, please read LICENCE file.


# --------------------------------
cmake_minimum_required(VERSION 2.6)

ENABLE_LANGUAGE(CXX)

FIND_PACKAGE( OpenMP REQUIRED )

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# --------------------------------
set(TRIMESH_PATH "${CMAKE_CURRENT_SOURCE_DIR}/trimesh2" CACHE PATH "Path to Trimesh library")
set(TRIMESH_LIB "${TRIMESH_PATH}/lib.Darwin64/libtrimesh.a" CACHE FILEPATH "Trimesh library")

message( STATUS "TRIMESH_PATH: " ${TRIMESH_PATH} )
message( STATUS "TRIMESH_LIB: " ${TRIMESH_LIB} )

# --------------------------------
add_definitions(-DUSE_SOS)

set(SOS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Detri_2.6.a" CACHE PATH "Path to SOS library")
set(SOS_LIB "${SOS_PATH}/build/lib/libSoS.a" CACHE FILEPATH "SOS library")

message( STATUS "SOS_PATH: " ${SOS_PATH} )
message( STATUS "SOS_LIB: " ${SOS_LIB} )

# --------------------------------

include_directories(
        ./include
        ${TRIMESH_PATH}/include
        ${SOS_PATH}/basic
        ${SOS_PATH}/lia
        ${SOS_PATH}/sos
)

set(SOURCE ./src/TriMeshJ.cpp ./src/JacobiSet.cpp ./src/main.cpp)
set(HEADERS ./include/TriMeshJ.h ./include/JacobiSet.h ./include/sos_utils.h)

add_executable(JacobiSetComputation ${SOURCE} ${HEADER})
target_link_libraries(JacobiSetComputation ${TRIMESH_LIB} ${SOS_LIB})
