cmake_minimum_required(VERSION 3.16)

project(LightFEM VERSION 0.6.0 LANGUAGES CXX)

option(LIGHTFEM_BUILD_TESTS "Build LightFEM tests" ON)

add_library(lightfem_core
    src/model.cpp
    src/input_parser.cpp
    src/validation.cpp
    src/element.cpp
    src/sparse.cpp
    src/reordering.cpp
    src/solver.cpp
    src/vtk_writer.cpp
)

target_include_directories(lightfem_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_features(lightfem_core PUBLIC cxx_std_17)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(lightfem_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(lightfem_core PUBLIC LIGHTFEM_HAS_OPENMP=1)
else()
    message(WARNING "OpenMP was not found: serial assembly and serial SpMV remain available")
endif()

if(MSVC)
    target_compile_options(lightfem_core PRIVATE /W4 /permissive- /utf-8)
else()
    target_compile_options(lightfem_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(lightfem src/main.cpp)
target_link_libraries(lightfem PRIVATE lightfem_core)

if(MSVC)
    target_compile_options(lightfem PRIVATE /W4 /permissive- /utf-8)
else()
    target_compile_options(lightfem PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(GNUInstallDirs)
install(TARGETS lightfem RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(LIGHTFEM_BUILD_TESTS)
    enable_testing()
    add_executable(lightfem_tests tests/test_input_validation.cpp)
    target_link_libraries(lightfem_tests PRIVATE lightfem_core)
    if(MSVC)
        target_compile_options(lightfem_tests PRIVATE /W4 /permissive- /utf-8)
    else()
        target_compile_options(lightfem_tests PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    add_test(NAME lightfem_input_validation COMMAND lightfem_tests)

    add_executable(lightfem_solver_tests tests/test_solver.cpp)
    target_link_libraries(lightfem_solver_tests PRIVATE lightfem_core)
    if(MSVC)
        target_compile_options(lightfem_solver_tests PRIVATE /W4 /permissive- /utf-8)
    else()
        target_compile_options(lightfem_solver_tests PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    add_test(NAME lightfem_solver COMMAND lightfem_solver_tests)
endif()
