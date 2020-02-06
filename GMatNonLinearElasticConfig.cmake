# GMatNonLinearElastic cmake module
#
# This module sets the target:
#
#     GMatNonLinearElastic
#
# In addition, it sets the following variables:
#
#     GMatNonLinearElastic_FOUND - true if GMatNonLinearElastic found
#     GMatNonLinearElastic_VERSION - GMatNonLinearElastic's version
#     GMatNonLinearElastic_INCLUDE_DIRS - the directory containing GMatNonLinearElastic headers
#
# The following support targets are defined to simplify things:
#
#     GMatNonLinearElastic::compiler_warnings - enable compiler warnings
#     GMatNonLinearElastic::assert - enable GMatNonLinearElastic assertions
#     GMatNonLinearElastic::debug - enable all assertions (slow)

include(CMakeFindDependencyMacro)

# Define target "GMatNonLinearElastic"

if(NOT TARGET GMatNonLinearElastic)
    include("${CMAKE_CURRENT_LIST_DIR}/GMatNonLinearElasticTargets.cmake")
    get_target_property(
        GMatNonLinearElastic_INCLUDE_DIRS
        GMatNonLinearElastic
        INTERFACE_INCLUDE_DIRECTORIES)
endif()

# Find dependencies

find_dependency(xtensor)

# Define support target "GMatNonLinearElastic::compiler_warnings"

if(NOT TARGET GMatNonLinearElastic::compiler_warnings)
    add_library(GMatNonLinearElastic::compiler_warnings INTERFACE IMPORTED)
    if(MSVC)
        set_property(
            TARGET GMatNonLinearElastic::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            /W4)
    else()
        set_property(
            TARGET GMatNonLinearElastic::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            -Wall -Wextra -pedantic -Wno-unknown-pragmas)
    endif()
endif()

# Define support target "GMatNonLinearElastic::assert"

if(NOT TARGET GMatNonLinearElastic::assert)
    add_library(GMatNonLinearElastic::assert INTERFACE IMPORTED)
    set_property(
        TARGET GMatNonLinearElastic::assert
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        GMATNONLINEARELASTIC_ENABLE_ASSERT)
endif()

# Define support target "GMatNonLinearElastic::debug"

if(NOT TARGET GMatNonLinearElastic::debug)
    add_library(GMatNonLinearElastic::debug INTERFACE IMPORTED)
    set_property(
        TARGET GMatNonLinearElastic::debug
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        XTENSOR_ENABLE_ASSERT GMATNONLINEARELASTIC_ENABLE_ASSERT)
endif()
