cmake_minimum_required(VERSION 3.10)
project(polygonRSA)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

find_package(OpenMP REQUIRED)

file(GLOB_RECURSE src CONFIGURE_DEPENDS "polygonRSA/*.cpp" "polygonRSA/*.h" "polygonRSA/*.tpp")
add_executable(polygonRSA ${src})
target_link_libraries(polygonRSA PRIVATE OpenMP::OpenMP_CXX)
