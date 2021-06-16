set(TARGET mdk)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake" ${CMAKE_MODULE_PATH})

add_library(${TARGET})

file(GLOB_RECURSE PUBLIC_SOURCES include/*)
file(GLOB_RECURSE PRIVATE_SOURCES src/*)

target_sources(${TARGET}
    PUBLIC ${PUBLIC_SOURCES}
    PRIVATE ${PRIVATE_SOURCES})

target_include_directories(${TARGET}
    PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/mdk>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/mdk>
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

target_compile_options(${TARGET}
    PRIVATE
    -Wall -Wextra -fopenmp -Wno-unused-parameter
    $<$<CONFIG:Debug>:-pg -Og>
    $<$<CONFIG:Release>:-Ofast>)

target_compile_features(${TARGET}
    PUBLIC cxx_std_17)

find_package(Eigen3 REQUIRED)
find_package(Catch2 REQUIRED)
find_package(OpenMP REQUIRED)

target_link_libraries(${TARGET}
    PUBLIC Eigen3::Eigen
    PRIVATE Catch2::Catch2
    PRIVATE OpenMP::OpenMP_CXX)

add_subdirectory(docs)