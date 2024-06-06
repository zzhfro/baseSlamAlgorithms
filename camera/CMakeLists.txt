cmake_minimum_required(VERSION 2.8)
project(camera)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-O3")
# 查找 Eigen 库

include_directories("/usr/include/eigen3")


# 添加可执行文件
add_executable(testcamera testcamera.cpp)




