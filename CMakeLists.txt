project(BasicSPH)
cmake_minimum_required(VERSION 2.8)

set(SRC ./src)
file(GLOB SOURCEFILES ${SRC}/*.cpp ${SRC}/*.h ${SRC}/*.hpp)

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SOURCEFILES})

