
# Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

# Set Compiler
set(CMAKE_C_COMPILER  "icl")
set(CMAKE_CXX_COMPILER  "icl")

# Project's name
project(SVMPath CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# --- OpenMP ---
include (FindOpenMP)
if (OPENMP_FOUND)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${OpenMP_LINKER_FLAGS}") 
endif()

# Blas
include( FindBLAS )
if (BLAS_FOUND)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BLAS_CXX_FLAGS}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${BLAS_C_FLAGS}")
  set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}") 
endif()

# Lapack
include( FindLAPACK )
if (FindLAPACK_FOUND)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LAPACK_CXX_FLAGS}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${LAPACK_C_FLAGS}")
  set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}") 
endif()

# Armadillo
set(ARMADILLO_INCLUDE_DIR "../Armadillo/include")
if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(ARMADILLO_LIBRARY "../Armadillo/lib/armadillod.lib")
else()
	set(ARMADILLO_LIBRARY "../Armadillo/lib/armadillo.lib")
endif()
include( FindArmadillo ) 

# Source Files
SET(TARGET_SOURCE "svm.cpp")

include_directories(${PROJECT_SOURCE_DIR} ${ARMADILLO_INCLUDE_DIRS})
add_executable(svmpath ${TARGET_SOURCE})
target_link_libraries(svmpath ${ARMADILLO_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
