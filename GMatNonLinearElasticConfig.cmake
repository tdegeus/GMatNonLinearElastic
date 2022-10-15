# GMatNonLinearElastic cmake module
#
# This module sets the target:
#
#     GMatNonLinearElastic
#
# In addition, it sets the following variables:
#
#     GMatNonLinearElastic_FOUND - true if the library is found
#     GMatNonLinearElastic_VERSION - the library's version
#     GMatNonLinearElastic_INCLUDE_DIRS - directory containing the library's headers
#
# The following support targets are defined to simplify things:
#
#     GMatNonLinearElastic::compiler_warnings - enable compiler warnings
#     GMatNonLinearElastic::assert - enable library assertions
#     GMatNonLinearElastic::debug - enable all assertions (slow)

include(CMakeFindDependencyMacro)

# Define target "GMatNonLinearElastic"

if(NOT TARGET GMatNonLinearElastic)
    include("${CMAKE_CURRENT_LIST_DIR}/GMatNonLinearElasticTargets.cmake")
endif()

# Define "GMatNonLinearElastic_INCLUDE_DIRS"

get_target_property(
    GMatNonLinearElastic_INCLUDE_DIRS
    GMatNonLinearElastic
    INTERFACE_INCLUDE_DIRECTORIES)

# Find dependencies

find_dependency(GMatElastic)
find_dependency(GMatTensor)
find_dependency(xtensor)

# Define support target "GMatNonLinearElastic::compiler_warnings"

if(NOT TARGET GMatNonLinearElastic::compiler_warnings)
    add_library(GMatNonLinearElastic::compiler_warnings INTERFACE IMPORTED)
    target_link_libraries(GMatNonLinearElastic::compiler_warnings INTERFACE
        GMatTensor::compiler_warnings)
endif()

# Define support target "GMatNonLinearElastic::assert"

if(NOT TARGET GMatNonLinearElastic::assert)
    add_library(GMatNonLinearElastic::assert INTERFACE IMPORTED)
    set_property(
        TARGET GMatNonLinearElastic::assert
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        GMATELASTIC_ENABLE_ASSERT
        GMATNONLINEARELASTIC_ENABLE_ASSERT)
endif()

# Define support target "GMatNonLinearElastic::debug"

if(NOT TARGET GMatNonLinearElastic::debug)
    add_library(GMatNonLinearElastic::debug INTERFACE IMPORTED)
    set_property(
        TARGET GMatNonLinearElastic::debug
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        GMATELASTIC_ENABLE_ASSERT
        GMATNONLINEARELASTIC_ENABLE_ASSERT
        GMATTENSOR_ENABLE_ASSERT
        XTENSOR_ENABLE_ASSERT)
endif()
