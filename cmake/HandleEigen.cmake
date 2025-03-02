###############################################################################
# Option for using system Eigen or GTSAM-bundled Eigen
option(GTSAM_USE_SYSTEM_EIGEN "Find and use system-installed Eigen. If 'off', use the one bundled with GTSAM" ON)

if(NOT GTSAM_USE_SYSTEM_EIGEN)
  # This option only makes sense if using the embedded copy of Eigen, it is
  # used to decide whether to *install* the "unsupported" module:
  option(GTSAM_WITH_EIGEN_UNSUPPORTED "Install Eigen's unsupported modules" OFF)
endif()

# Switch for using system Eigen or GTSAM-bundled Eigen
if(GTSAM_USE_SYSTEM_EIGEN)
    # Since Eigen 3.3.0 a Eigen3Config.cmake is available so use it.
    find_package(Eigen3 CONFIG REQUIRED) # need to find again as REQUIRED

    # The actual include directory (for BUILD cmake target interface):
    # Note: EIGEN3_INCLUDE_DIR points to some random location on some eigen
    # versions.  So here I use the target itself to get the proper include
    # directory (it is generated by cmake, thus has the correct path)
    get_target_property(GTSAM_EIGEN_INCLUDE_FOR_BUILD Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)

    # check if MKL is also enabled - can have one or the other, but not both!
    # Note: Eigen >= v3.2.5 includes our patches
    if(EIGEN_USE_MKL_ALL AND (EIGEN3_VERSION VERSION_LESS 3.2.5))
      message(FATAL_ERROR "MKL requires at least Eigen 3.2.5, and your system appears to have an older version. Disable GTSAM_USE_SYSTEM_EIGEN to use GTSAM's copy of Eigen, or disable GTSAM_WITH_EIGEN_MKL")
    endif()

    # Check for Eigen version which doesn't work with MKL
    # See http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1527 for details.
    if(EIGEN_USE_MKL_ALL AND (EIGEN3_VERSION VERSION_EQUAL 3.3.4))
        message(FATAL_ERROR "MKL does not work with Eigen 3.3.4 because of a bug in Eigen. See http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1527. Disable GTSAM_USE_SYSTEM_EIGEN to use GTSAM's copy of Eigen, disable GTSAM_WITH_EIGEN_MKL, or upgrade/patch your installation of Eigen.")
    endif()

    message("USING EIGEN SYSTEM INSTALLATION")
else()
    # Use bundled Eigen include path.
    # Clear any variables set by FindEigen3
    if(EIGEN3_INCLUDE_DIR)
        set(EIGEN3_INCLUDE_DIR NOTFOUND CACHE STRING "" FORCE)
    endif()

    # set full path to be used by external projects
    # this will be added to GTSAM_INCLUDE_DIR by gtsam_extra.cmake.in
    set(GTSAM_EIGEN_INCLUDE_FOR_INSTALL "include/gtsam/3rdparty/Eigen/")

    # The actual include directory (for BUILD cmake target interface):
    set(GTSAM_EIGEN_INCLUDE_FOR_BUILD "${GTSAM_SOURCE_DIR}/gtsam/3rdparty/Eigen")

    add_library(gtsam_eigen3 INTERFACE)

    target_include_directories(gtsam_eigen3 INTERFACE
      $<BUILD_INTERFACE:${GTSAM_EIGEN_INCLUDE_FOR_BUILD}>
      $<INSTALL_INTERFACE:${GTSAM_EIGEN_INCLUDE_FOR_INSTALL}>
    )
    add_library(Eigen3::Eigen ALIAS gtsam_eigen3)

    install(TARGETS gtsam_eigen3 EXPORT GTSAM-exports PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    list(APPEND GTSAM_EXPORTED_TARGETS gtsam_eigen3)
    set(GTSAM_EXPORTED_TARGETS "${GTSAM_EXPORTED_TARGETS}")
    message("USING EIGEN BUNDLED IN GTSAM")
endif()

# Detect Eigen version:
set(EIGEN_VER_H "${GTSAM_EIGEN_INCLUDE_FOR_BUILD}/Eigen/src/Core/util/Macros.h")
if (EXISTS ${EIGEN_VER_H})
    file(READ "${EIGEN_VER_H}" STR_EIGEN_VERSION)

    # Extract the Eigen version from the Macros.h file, lines "#define EIGEN_WORLD_VERSION  XX", etc...

    string(REGEX MATCH "EIGEN_WORLD_VERSION[ ]+[0-9]+" GTSAM_EIGEN_VERSION_WORLD "${STR_EIGEN_VERSION}")
    string(REGEX MATCH "[0-9]+" GTSAM_EIGEN_VERSION_WORLD "${GTSAM_EIGEN_VERSION_WORLD}")

    string(REGEX MATCH "EIGEN_MAJOR_VERSION[ ]+[0-9]+" GTSAM_EIGEN_VERSION_MAJOR "${STR_EIGEN_VERSION}")
    string(REGEX MATCH "[0-9]+" GTSAM_EIGEN_VERSION_MAJOR "${GTSAM_EIGEN_VERSION_MAJOR}")

    string(REGEX MATCH "EIGEN_MINOR_VERSION[ ]+[0-9]+" GTSAM_EIGEN_VERSION_MINOR "${STR_EIGEN_VERSION}")
    string(REGEX MATCH "[0-9]+" GTSAM_EIGEN_VERSION_MINOR "${GTSAM_EIGEN_VERSION_MINOR}")

    set(GTSAM_EIGEN_VERSION "${GTSAM_EIGEN_VERSION_WORLD}.${GTSAM_EIGEN_VERSION_MAJOR}.${GTSAM_EIGEN_VERSION_MINOR}")

    message(STATUS "Found Eigen version: ${GTSAM_EIGEN_VERSION}")
else()
    message(WARNING "Cannot determine Eigen version, missing file: `${EIGEN_VER_H}`")
endif ()

if (MSVC)
    if (GTSAM_SHARED_LIB)
        # mute eigen static assert to avoid errors in shared lib
        list_append_cache(GTSAM_COMPILE_DEFINITIONS_PUBLIC EIGEN_NO_STATIC_ASSERT)
    endif()
    list_append_cache(GTSAM_COMPILE_OPTIONS_PRIVATE "/wd4244") # Disable loss of precision which is thrown all over our Eigen
endif()
