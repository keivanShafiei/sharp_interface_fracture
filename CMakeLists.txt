cmake_minimum_required(VERSION 3.14)
project(SharpInterfaceFracture CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ============================================================================
# Compiler flags
# ============================================================================
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -Wall -Wextra")

# ============================================================================
# Find dependencies
# ============================================================================

# PETSc (required for linear system solver)
find_package(PETSc REQUIRED)
# If PETSc not found via cmake module, try pkg-config:
if(NOT PETSc_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(PETSC REQUIRED petsc)
    set(PETSC_INCLUDES ${PETSC_INCLUDE_DIRS})
    set(PETSC_LIBRARIES ${PETSC_LINK_LIBRARIES})
endif()

# Eigen3 (header-only, required for matrix operations)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# CGAL (required for constrained Delaunay remeshing)
find_package(CGAL REQUIRED COMPONENTS Core)

# ============================================================================
# Main executable
# ============================================================================
set(SOURCES
    src/main.cpp
    src/assembler.cpp
    src/config_force.cpp
    src/mesh.cpp
    src/crack.cpp
    src/field_transfer.cpp
    src/remesher.cpp
    src/vtk.cpp
    src/simulation_config.cpp
)

add_executable(sharp_fracture ${SOURCES})

target_include_directories(sharp_fracture PRIVATE
    include
    ${PETSC_INCLUDES}
)

target_link_libraries(sharp_fracture PRIVATE
    ${PETSC_LIBRARIES}
    Eigen3::Eigen
    CGAL::CGAL
    CGAL::CGAL_Core
)

# ============================================================================
# Tests
# ============================================================================
enable_testing()

# Test 1: Patch test
add_executable(test_patch_test
    tests/test_patch_test.cpp
    src/assembler.cpp
    src/mesh.cpp
    src/crack.cpp
    src/config_force.cpp
    src/field_transfer.cpp
    src/remesher.cpp
    src/vtk.cpp
    src/simulation_config.cpp
)
target_include_directories(test_patch_test PRIVATE include ${PETSC_INCLUDES})
target_link_libraries(test_patch_test PRIVATE
    ${PETSC_LIBRARIES} Eigen3::Eigen CGAL::CGAL CGAL::CGAL_Core)
add_test(NAME PatchTest COMMAND test_patch_test)

# Test 2: Elasticity tensor consistency
add_executable(test_elasticity
    tests/test_elasticity_consistency.cpp
    src/config_force.cpp
    src/mesh.cpp
)
target_include_directories(test_elasticity PRIVATE include ${PETSC_INCLUDES})
target_link_libraries(test_elasticity PRIVATE
    ${PETSC_LIBRARIES} Eigen3::Eigen CGAL::CGAL CGAL::CGAL_Core)
add_test(NAME ElasticityConsistency COMMAND test_elasticity)

# ============================================================================
# Installation
# ============================================================================
install(TARGETS sharp_fracture DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)
