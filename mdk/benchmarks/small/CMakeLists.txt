set(TARGET bench-small)
add_executable(${TARGET} main.cpp)

target_link_libraries(${TARGET}
    PRIVATE mdk)

target_compile_features(${TARGET}
    PRIVATE cxx_std_17)
