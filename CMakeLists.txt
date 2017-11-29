
# Define Compiler
set(CMAKE_C_COMPILER  "icl")
set(CMAKE_CXX_COMPILER  "icl")

# Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

# Project's name
project(SVMPath LANGUAGES C CXX)

# --- OpenMP ---
if(NOT HAVE_TBB AND NOT HAVE_CSTRIPES)
  include (FindOpenMP)
   if (OPENMP_FOUND)
      set (HAVE_OPENMP TRUE)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   else()
      set ( HAVE_OPENMP FALSE)
   endif()
else()
set(HAVE_OPENMP 0)
endif()

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

SET(TARGET_SOURCE "svm.cpp")

# 
add_executable(svm-train svm-train.c ${TARGET_SOURCE})

add_executable(svm-scale svm-scale.c ${TARGET_SOURCE})

add_executable(svm-predict svm-predict.c ${TARGET_SOURCE})


##
