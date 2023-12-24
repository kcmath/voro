cmake_minimum_required(VERSION 3.20)
project(VoroMesh VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)

find_package (VTK REQUIRED PATHS "")
find_package (Eigen3 REQUIRED PATHS "")
set(VORO_INCLUDE_DIR "")
set(VORO_LIBRARY "")

set(EIGEN_INCLUDE_DIR "")


file(GLOB_RECURSE MESH_SOURCES "src/*.cpp" "src/*/*.cpp")
add_library(VoroMesh STATIC ${MESH_SOURCES})

target_include_directories(VoroMesh PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    ${VORO_INCLUDE_DIR}
    $<INSTALL_INTERFACE:include>
    ${EIGEN_INCLUDE_DIR}
    )
    
link_directories("")    
target_link_libraries (VoroMesh ${VTK_LIBRARIES} ${VORO_LIBRARY})

add_executable(driver driver/driver.cpp)
target_link_libraries(driver PUBLIC VoroMesh)
