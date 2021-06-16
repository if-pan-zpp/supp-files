set(TARGET 1ubq)
add_executable(${TARGET} main.cpp)

target_link_libraries(${TARGET}
    PRIVATE mdk)

target_compile_features(${TARGET}
    PRIVATE cxx_std_17)

add_custom_command(TARGET ${TARGET} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_LIST_DIR}/data $<TARGET_FILE_DIR:${TARGET}>/data)