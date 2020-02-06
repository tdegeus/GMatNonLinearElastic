cmake_minimum_required(VERSION 3.0)

# Basic settings
# ==============

project(GMatNonLinearElastic)

option(BUILD_TESTS "${PROJECT_NAME} Build tests" OFF)

# Version
# =======

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/GMatNonLinearElastic/config.h" GMatNonLinearElastic_version_defines
     REGEX "#define GMATNONLINEARELASTIC_VERSION_(MAJOR|MINOR|PATCH)")

foreach(ver ${GMatNonLinearElastic_version_defines})
    if(ver MATCHES "#define GMATNONLINEARELASTIC_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
    set(GMATNONLINEARELASTIC_VERSION_${CMAKE_MATCH_1}
        "${CMAKE_MATCH_2}"
        CACHE INTERNAL "")
    endif()
endforeach()

set(GMATNONLINEARELASTIC_VERSION
    ${GMATNONLINEARELASTIC_VERSION_MAJOR}.${GMATNONLINEARELASTIC_VERSION_MINOR}.${GMATNONLINEARELASTIC_VERSION_PATCH})

message(STATUS "Building GMatNonLinearElastic v${GMATNONLINEARELASTIC_VERSION}")

# Set target
# ==========

find_package(xtensor REQUIRED)

add_library(GMatNonLinearElastic INTERFACE)

target_include_directories(GMatNonLinearElastic INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(GMatNonLinearElastic INTERFACE xtensor)

# Installation
# ============

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION include)

install(TARGETS GMatNonLinearElastic EXPORT GMatNonLinearElastic-targets)

install(
  EXPORT GMatNonLinearElastic-targets
  FILE GMatNonLinearElasticTargets.cmake
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GMatNonLinearElastic")

set(_GMATNONLINEARELASTIC ${CMAKE_SIZEOF_VOID_P})
unset(CMAKE_SIZEOF_VOID_P)

write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/GMatNonLinearElasticConfigVersion.cmake"
    VERSION ${GMATNONLINEARELASTIC_VERSION}
    COMPATIBILITY AnyNewerVersion)

set(CMAKE_SIZEOF_VOID_P ${_GMATNONLINEARELASTIC})

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/GMatNonLinearElasticConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/GMatNonLinearElasticConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GMatNonLinearElastic")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GMatNonLinearElastic.pc.in"
               "${CMAKE_CURRENT_BINARY_DIR}/GMatNonLinearElastic.pc" @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/GMatNonLinearElastic.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")

# Add builds
# ==========

include("GMatNonLinearElasticConfig.cmake")

if(BUILD_TESTS)
    add_subdirectory(test)
endif()
