# Python and MATLAB Wrapper Guide

The optional wrappers use `src/wrap_interface.i`. The generated top namespace is `estimation_gears`.

The pinned `lib/wrap` revision `fc811c1` supports the Python wrapper but predates the MATLAB fixes for constant string references, fixed-width integers, and MEX error recovery. Run the MATLAB regression matrix with an explicit gtwrap checkout containing merge `1c27f1f` or later. Moving the repository gitlink remains a separate maintenance operation.

## Options

- `EstimationGears_for_SpaceNav_BUILD_PYTHON_WRAPPER`
- `EstimationGears_for_SpaceNav_BUILD_MATLAB_WRAPPER`
- `EstimationGears_for_SpaceNav_GTWRAP_ROOT_DIR`
- `EstimationGears_for_SpaceNav_WRAPPER_INTERFACE_FILES`
- `EstimationGears_for_SpaceNav_GTWRAP_TOP_NAMESPACE`
- `EstimationGears_for_SpaceNav_GTWRAP_RUNTIME_DEPENDENCY_TARGETS`

`build_lib.sh -p` and `-m` set the corresponding defaults. Wrapper checkout updates and submodule initialization remain off unless explicit maintenance flags are provided.

## Python

```bash
./build_lib.sh -N -p --gtwrap-root lib/wrap
ctest --test-dir build --output-on-failure -L python
cd python
python -m pip wheel . --no-build-isolation --no-deps
```

CMake configures an immutable build-tree copy of the source package. Wheel creation copies only the declared extension module and exact project-owned shared runtime artifacts. Loader-relative runtime paths allow the extension to find those colocated libraries.

`python/EstimationGears_for_SpaceNav/_wrapper_build.py` is checkout-only metadata. It must not be installed or included in wheels.

## MATLAB

```bash
./build_lib.sh -N -m --gtwrap-root /path/to/fixed/wrap
ctest --test-dir build --output-on-failure -L matlab
```

The regression cases cover object lifetime, strings by value and constant reference, unsigned identifiers, error recovery, standard output, `clear all`, and optional tcmalloc linkage.

CMake Python and MATLAB installation destinations stay relative to `CMAKE_INSTALL_PREFIX`; pip owns installation into active Python environments.
