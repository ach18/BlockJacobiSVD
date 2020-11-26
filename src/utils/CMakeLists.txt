cmake_minimum_required(VERSION 3.12)
set(SOURCES
    util.cpp
    matrix.cpp
    block.cpp
    )

# Check if Boost is already installed
find_package(Boost 1.70)

if(NOT Boost_FOUND)
    # If no installation is found add it
    add_subdirectory(boost-cmake)
else()
    # If an installation is found, use it
    set(BOOST_SOURCE ${Boost_INCLUDE_DIRS})
endif()

add_library(utils-lib ${SOURCES})
target_include_directories(utils-lib PUBLIC .)
message("boost source ${BOOST_SOURCE}")
target_include_directories(utils-lib PUBLIC ${BOOST_SOURCE})
