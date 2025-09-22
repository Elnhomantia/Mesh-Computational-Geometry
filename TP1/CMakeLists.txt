cmake_minimum_required(VERSION 3.16)

project(TP1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)

set(INCLUDES_FILES
    include/point.h
    include/face.h
    include/vertex.h
    include/mesh.h
    include/circulatorOnFaces.h
    include/circulatorOnVertices.h
)

set(SOURCES_FILES
    src/OFFparser.cpp
    src/meshChecker.cpp
    src/mesh.cpp
    src/circulatorOnFaces.cpp
    src/circulatorOnVertices.cpp
    main.cpp
)

add_executable(TP1 ${INCLUDES_FILES} ${SOURCES_FILES})

target_include_directories(TP1 PUBLIC include)

set(RESOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
file(GLOB RESOURCE_FILES "${RESOURCE_DIR}/*.off")
foreach(f ${RESOURCE_FILES})
    configure_file(${f} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
endforeach()