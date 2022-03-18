cmake_minimum_required(VERSION 3.17)
project(fractals)

set(CMAKE_CXX_STANDARD 14)


include(ExternalProject)

ExternalProject_Add(
    eigen
    PREFIX "lib/eigen"
    GIT_REPOSITORY "https://gitlab.com/libeigen/eigen"
    TIMEOUT 10
    CMAKE_ARGS ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    UPDATE_COMMAND ""
)



ExternalProject_Get_Property(eigen source_dir)
set(EIGEN_INCLUDE_DIR ${source_dir}/)
set_target_properties(eigen PROPERTIES EXCLUDE_FROM_ALL TRUE)
include_directories(
    ${EIGEN_INCLUDE_DIR}
)
add_executable(fractals mandelbrot.cpp)
add_dependencies(fractals eigen)