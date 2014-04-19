include_directories(${BACTELIZE_SOURCE_DIR}/src)
link_directories(${BACTELIZE_BINARY_DIR}/src)

ADD_EXECUTABLE(bactelize bactelize.cxx)
TARGET_LINK_LIBRARIES(bactelize ${ITK_LIBRARIES})

