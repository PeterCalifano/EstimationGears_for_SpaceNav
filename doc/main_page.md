# EstimationGears_for_SpaceNav {#mainpage}

EstimationGears provides MATLAB estimation algorithms and a C++20/CUDA library scaffold for spacecraft-navigation applications.

## Entry points

- [Repository overview](../README.md)
- [C++ and CUDA builds](cpp_cuda_build.md)
- [Build script](build_script_doc.md)
- [Python and MATLAB wrappers](wrappers.md)
- [Testing and CI](testing_and_ci.md)
- [Versioning](versioning.md)
- [Logging](logging.md)
- [Documentation workflow](documentation_workflow.md)

## Installed CMake API

```cmake
find_package(EstimationGears_for_SpaceNav REQUIRED CONFIG)
target_link_libraries(
    application
    PRIVATE EstimationGears_for_SpaceNav::EstimationGears_for_SpaceNav)
```

The native build is independent of MATLAB-only source/runtime checkouts under `lib/`. CUDA, Python wrapping, and MATLAB wrapping are optional.
