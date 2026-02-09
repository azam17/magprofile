cmake_minimum_required(VERSION 3.16)
project(HalalSeq C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ---- Core library (everything except main.c) -----------------------
file(GLOB CORE_SRCS ${CMAKE_SOURCE_DIR}/src/*.c)
list(REMOVE_ITEM CORE_SRCS ${CMAKE_SOURCE_DIR}/src/main.c)

add_library(halalseq_core STATIC ${CORE_SRCS})
target_include_directories(halalseq_core PUBLIC
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib)
target_compile_options(halalseq_core PRIVATE -Wall -Wextra)
target_link_libraries(halalseq_core PRIVATE z m)

# ---- CLI binary ----------------------------------------------------
add_executable(halalseq ${CMAKE_SOURCE_DIR}/src/main.c)
target_link_libraries(halalseq PRIVATE halalseq_core z m)

# ---- GUI binary (only if SDL2 is available) ------------------------
find_package(SDL2 QUIET)
if(SDL2_FOUND)
    add_executable(HalalSeqGUI WIN32 MACOSX_BUNDLE
        ${CMAKE_SOURCE_DIR}/gui/main_gui.c
        ${CMAKE_SOURCE_DIR}/gui/gui_analysis.c
        ${CMAKE_SOURCE_DIR}/lib/tinyfiledialogs.c)

    target_include_directories(HalalSeqGUI PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/lib
        ${CMAKE_SOURCE_DIR}/gui)

    target_link_libraries(HalalSeqGUI PRIVATE
        halalseq_core
        SDL2::SDL2
        z m)

    if(APPLE)
        target_link_libraries(HalalSeqGUI PRIVATE "-framework Cocoa")
        set_target_properties(HalalSeqGUI PROPERTIES
            MACOSX_BUNDLE_INFO_PLIST ${CMAKE_SOURCE_DIR}/gui/Info.plist
            MACOSX_BUNDLE_BUNDLE_NAME "HalalSeq"
            MACOSX_BUNDLE_BUNDLE_VERSION "0.1.0"
            MACOSX_BUNDLE_SHORT_VERSION_STRING "0.1.0")
    endif()

    if(WIN32)
        target_link_libraries(HalalSeqGUI PRIVATE comdlg32 ole32)
    endif()

    message(STATUS "SDL2 found — GUI target enabled")
else()
    message(STATUS "SDL2 not found — GUI target disabled (CLI only)")
endif()
