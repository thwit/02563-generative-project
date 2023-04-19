cmake_minimum_required(VERSION 3.26)

project(erosionproject)

# Find the Eigen3 library
find_package(Eigen3 CONFIG REQ  UIRED)

# Create an executable target
add_executable(main erosion.cpp)

# Link the Eigen3 library to the main executable
target_link_libraries(main PRIVATE Eigen3::Eigen)