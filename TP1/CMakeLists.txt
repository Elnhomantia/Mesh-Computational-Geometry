cmake_minimum_required(VERSION 3.16)

project(TP1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(INCLUDES_FILES
    include/point.h
    include/triangle.h
    include/vertex.h
)

set(SOURCES_FILES
    OFFparser.cpp    
    main.cpp
)

add_executable(TP1 ${INCLUDES_FILES} ${SOURCES_FILES})

target_include_directories(TP1 PUBLIC include)

set(RESOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
file(GLOB RESOURCE_FILES "${RESOURCE_DIR}/*.off")
foreach(f ${RESOURCE_FILES})
    configure_file(${f} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
endforeach()