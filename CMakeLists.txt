cmake_minimum_required(VERSION 3.20)

project(cmake_test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(reactor
    src/main.cpp
    # More .cpp if necessary    
)

target_include_directories(reactor
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/headers
        ${CMAKE_CURRENT_SOURCE_DIR}/source-terms
)