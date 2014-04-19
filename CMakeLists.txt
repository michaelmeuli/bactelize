cmake_minimum_required(VERSION 2.6)
 
PROJECT(BACTELIZE)
 
FIND_PACKAGE(ITK REQUIRED SCIFIO)
INCLUDE(${ITK_USE_FILE})
 
if(NOT SCIFIO_LOADED)
  message(FATAL_ERROR "SCIFIO module is required but not available.")
endif()
 
add_subdirectory(src)

