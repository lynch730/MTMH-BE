cmake_minimum_required(VERSION ${CMAKE_MINIMUM_REQUIRED_VERSION})


get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})
project(${ProjectId})




include_directories(${ProjectId} ${__INCLUDE_DIRECTORIES})
add_executable(${ProjectId} ${ProjectId}.cpp)
target_link_libraries(${ProjectId} ${__LINK_LIBRARIES})
