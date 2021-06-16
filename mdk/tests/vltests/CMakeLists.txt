add_executable(vltests main.cpp)

target_compile_options(vltests
    PRIVATE
    -march=native -fopenmp
    $<$<CONFIG:Debug>:-Og -g>
    $<$<CONFIG:Release>:-Ofast>)

find_package(OpenMP)
target_link_libraries(vltests
    PRIVATE OpenMP::OpenMP_CXX)