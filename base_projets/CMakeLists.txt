# CMake entry point
cmake_minimum_required (VERSION 2.6)
project(Base_Projet)


#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
#et(CMAKE_CXX_COMPILER "g++-7")

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)

find_package(OpenGL REQUIRED)

add_subdirectory (./external)

include_directories(src/)

include_directories(
        external/glfw-2.7.6/include/
        external/glm-0.9.4.0/
        external/glew-1.9.0/include/
        external/eigen-3.3.4/
)


set(ALL_LIBS
	GLFW_276
        GLEW_190
)

add_definitions(
	-DTW_STATIC
	-DTW_NO_LIB_PRAGMA
	-DTW_NO_DIRECT3D
	-DGLEW_STATIC
	-D_CRT_SECURE_NO_WARNINGS

)



file(
    GLOB
    shader_file
    shader/*.glsl
    )

file(
    GLOB
    source_files
    src/AbstractMesh.cpp
    src/DataSet.cpp
    src/Graph.cpp
    src/ImplicitFunction.cpp
    src/Mesh.cpp
    src/MeshHE.cpp
    src/Object.cpp
    src/Plane.cpp
    src/shader.cpp
    src/*.h
)

add_executable(reconstruction src/main.cpp ${shader_file} ${source_files} )
add_executable(visu src/visuCloud.cpp ${shader_file} ${source_files})
target_link_libraries(visu ${ALL_LIBS})
target_link_libraries(reconstruction ${ALL_LIBS})
