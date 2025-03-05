cmake_minimum_required(VERSION 3.28)
project(UCSD_Renderer)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(TBB REQUIRED)

add_executable(UCSD_Renderer
        src/main.cpp
        src/util/VecMath.cpp
        src/Shape.cpp
        src/Scene.cpp
        src/Renderer.cpp)

target_include_directories(UCSD_Renderer PRIVATE C:/MyPC/Libs/FreeImage)

target_link_libraries(UCSD_Renderer PRIVATE C:/MyPC/libs/FreeImage/FreeImage.lib)
target_link_libraries(UCSD_Renderer PUBLIC TBB::tbb)
