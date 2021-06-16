cmake_minimum_required(VERSION 3.17)
project(mdk LANGUAGES CXX)

option(LEGACY_MODE "Work as close to f77 version as possible." OFF)
if (LEGACY_MODE)
  add_definitions(-DLEGACY_MODE)
endif()

add_subdirectory(mdk)
add_subdirectory(examples)
add_subdirectory(tests)
