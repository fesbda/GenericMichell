cmake_minimum_required(VERSION 3.18)

project(GenericMichell VERSION 1.0.0 LANGUAGES CXX)

include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(PYBIND11_FINDPYTHON ON)

option(MICHELL_USE_OPENMP "Enable OpenMP acceleration when available" ON)
option(MICHELL_BUILD_TESTS "Enable the Python smoke test" ON)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

FetchContent_Declare(
    pybind11
    URL https://github.com/pybind/pybind11/archive/refs/tags/v2.13.6.tar.gz
    URL_HASH SHA256=e08cb87f4773da97fa7b5f035de8763abc656d87d5773e62f6da0587d1f0ec20
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(pybind11)

FetchContent_Declare(
    glm
    URL https://github.com/g-truc/glm/archive/refs/tags/1.0.1.tar.gz
    URL_HASH SHA256=9f3174561fd26904b23f0db5e560971cbf9b3cbda0b280f04d5c379d03bf234c
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(glm)

pybind11_add_module(michell
    python/main.cpp
    ext/xfoil/XFoil.cpp
)

target_include_directories(michell PRIVATE src ext/xfoil)
target_compile_definitions(michell PRIVATE MICHELL_VERSION="${PROJECT_VERSION}")

if(TARGET glm::glm)
    target_link_libraries(michell PRIVATE glm::glm)
elseif(TARGET glm)
    target_link_libraries(michell PRIVATE glm)
else()
    message(FATAL_ERROR "GLM was configured, but no GLM target is available")
endif()

if(MICHELL_USE_OPENMP)
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(michell PRIVATE OpenMP::OpenMP_CXX)
    else()
        message(STATUS "OpenMP not found; building without OpenMP")
    endif()
endif()

target_compile_options(michell PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-fsized-deallocation>
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall>
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wextra>
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

install(TARGETS michell LIBRARY DESTINATION . RUNTIME DESTINATION .)

if(MICHELL_BUILD_TESTS)
    enable_testing()
    add_test(
        NAME smoke
        COMMAND ${CMAKE_COMMAND} -E env
                "PYTHONPATH=$<TARGET_FILE_DIR:michell>"
                ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke_test.py
    )
endif()
