cmake_minimum_required(VERSION "3.16.6")

#set(CMAKE_CXX_COMPILER g++)
#set(CMAEK_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a")

project(app)

#target_include_directories(C:/mtl4" "C:/boost_1_72_0")
include_directories("C:/mtl4" "C:/boost")

add_executable("${PROJECT_NAME}" poisson.cpp)



#set(CMAKE_C_COMPILER "C:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin")
#set(CMAKE_CXX_COMPILER "C:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin")