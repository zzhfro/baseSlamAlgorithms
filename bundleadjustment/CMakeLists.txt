cmake_minimum_required(VERSION 2.8)
project(sfm)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-O3")
# 查找 Eigen 库

include_directories("/usr/include/eigen3")
find_package(OpenCV REQUIRED)

# 添加可执行文件
#add_library(ep_shared SHARED epipolarAlgorithm.cpp)
add_executable(testba batest.cpp)
target_link_libraries(testba ${OpenCV_LIBS})
#target_link_libraries(testep ep_shared)