
cmake_minimum_required(VERSION 3.0)

# Options
# =======

project(GMatNonLinearElastic)

option(BUILD_TESTS "${PROJECT_NAME} Build tests" OFF)

# Internals
# =========

set(GMATNONLINEARELASTIC_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Version
# =======

file(
    STRINGS "${GMATNONLINEARELASTIC_INCLUDE_DIR}/GMatNonLinearElastic/config.h" GMatNonLinearElastic_version_defines
    REGEX "#define GMATNONLINEARELASTIC_VERSION_(MAJOR|MINOR|PATCH)")

foreach(ver ${GMatNonLinearElastic_version_defines})
    if(ver MATCHES "#define GMATNONLINEARELASTIC_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(GMATNONLINEARELASTIC_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
    endif()
endforeach()

set(GMATNONLINEARELASTIC_VERSION
    ${GMATNONLINEARELASTIC_VERSION_MAJOR}.${GMATNONLINEARELASTIC_VERSION_MINOR}.${GMATNONLINEARELASTIC_VERSION_PATCH})

message(STATUS "Building GMatNonLinearElastic v${GMATNONLINEARELASTIC_VERSION}")

# Build
# =====

set(GMATNONLINEARELASTIC_HEADERS
    ${GMATNONLINEARELASTIC_INCLUDE_DIR}/GMatNonLinearElastic/Cartesian3d.h
    ${GMATNONLINEARELASTIC_INCLUDE_DIR}/GMatNonLinearElastic/Cartesian3d.hpp
    ${GMATNONLINEARELASTIC_INCLUDE_DIR}/GMatNonLinearElastic/Cartesian3d_NonLinearElastic.hpp
    ${GMATNONLINEARELASTIC_INCLUDE_DIR}/GMatNonLinearElastic/Cartesian3d_Matrix.hpp
    ${GMATNONLINEARELASTIC_INCLUDE_DIR}/GMatNonLinearElastic/config.h)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(
    ${PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${GMATNONLINEARELASTIC_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>)

# Installation
# ============

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}-targets)

export(
    EXPORT ${PROJECT_NAME}-targets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")

install(
    FILES ${GMATNONLINEARELASTIC_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

set(GMATNONLINEARELASTIC_CMAKECONFIG_INSTALL_DIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE
    STRING "install path for ${PROJECT_NAME}Config.cmake")

configure_package_config_file(
    "GMatNonLinearElasticConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION ${GMATNONLINEARELASTIC_CMAKECONFIG_INSTALL_DIR})

set(_GMATNONLINEARELASTIC_CMAKE_SIZEOF_VOID_P ${CMAKE_SIZEOF_VOID_P})

unset(CMAKE_SIZEOF_VOID_P)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${GMATNONLINEARELASTIC_VERSION}
    COMPATIBILITY AnyNewerVersion)

set(CMAKE_SIZEOF_VOID_P ${_GMATNONLINEARELASTIC_CMAKE_SIZEOF_VOID_P})

install(
    FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION ${GMATNONLINEARELASTIC_CMAKECONFIG_INSTALL_DIR})

install(
    EXPORT ${PROJECT_NAME}-targets
    FILE "${PROJECT_NAME}Targets.cmake"
    DESTINATION ${GMATNONLINEARELASTIC_CMAKECONFIG_INSTALL_DIR})

configure_file(
    "GMatNonLinearElastic.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
     @ONLY)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")

# Add builds
# ==========

if(BUILD_TESTS)
    add_subdirectory(test)
endif()
