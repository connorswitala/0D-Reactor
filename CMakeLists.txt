cmake_minimum_required(VERSION 3.20)

project(cmake_test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(reactor
    src/5species.cpp
    # More .cpp if necessary    
)

add_executable(N2reactor
    src/2species.cpp
)

target_include_directories(reactor
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/includes/common
        ${CMAKE_CURRENT_SOURCE_DIR}/includes/physics
        ${CMAKE_CURRENT_SOURCE_DIR}/includes/xml
)