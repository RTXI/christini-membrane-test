cmake_minimum_required(VERSION 3.14)

project(
    membrane_test
    VERSION 0.1.0
    DESCRIPTION "measure pipette resistance and membrane properties during patch clamp experiments"
    HOMEPAGE_URL "https://rtxi.org/"
    LANGUAGES CXX
)

# These lines help with third-party tooling integration
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(RTXI_CMAKE_SCRIPTS "" CACHE PATH "In development mode, RTXI uses this to expose conan dependencies")
set(CMAKE_PREFIX_PATH "${RTXI_CMAKE_SCRIPTS}")
set(RTXI_PACKAGE_PATH "/usr/local/" CACHE PATH "Path hint to RTXI package information")
set(CMAKE_INSTALL_RPATH "${RTXI_PACKAGE_PATH}/lib")
list(APPEND CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})

# ---- find libraries ----
find_package(rtxi REQUIRED HINTS ${RTXI_PACKAGE_PATH})
find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets HINTS ${RTXI_CMAKE_SCRIPTS})
find_package(fmt REQUIRED)
find_package(GSL REQUIRED)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)


#################################################################################################
### You can modify within this region for finding additional libraries and linking them to your #
### custom plugin. Make sure to install them prior to configuration or else build will fail     #
### with linking and include errors!                                                            # 
#################################################################################################
add_library(
    membrane_test MODULE
    widget.cpp
    widget.hpp
)

# Consult library website for how to link them to your plugin using cmake
target_link_libraries(membrane_test PUBLIC 
    rtxi::rtxi rtxi::rtxidsp rtxi::rtxigen rtxi::rtxififo Qt5::Core Qt5::Gui Qt5::Widgets
    GSL::gsl
    dl fmt::fmt 
)

################################################################################################ 

# This will tell cmake to automatically scan .ui files and generate the corresponding header file
#set_property(TARGET membrane_test PROPERTY AUTOUIC ON)

# We need to tell cmake to use the c++ version used to compile the dependent library or else...
get_target_property(REQUIRED_COMPILE_FEATURE rtxi::rtxi INTERFACE_COMPILE_FEATURES)
target_compile_features(membrane_test PRIVATE ${REQUIRED_COMPILE_FEATURE})

install(
    TARGETS membrane_test
    DESTINATION ${RTXI_PACKAGE_PATH}/bin/rtxi_modules
)

