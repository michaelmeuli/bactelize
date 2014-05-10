cmake_minimum_required(VERSION 2.8)
 
project(bactelize)
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()

if(NOT SCIFIO_LOADED)
  message(FATAL_ERROR "SCIFIO module is required but not available.")
endif() 

add_executable(bactelize MACOSX_BUNDLE bactelize.cxx)
 
if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(bactelize ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(bactelize ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )

