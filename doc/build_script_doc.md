# `build_lib.sh` Guide

`build_lib.sh` is the supported native convenience entry point. It configures with CMake, builds through `cmake --build`, runs CTest by default, and optionally installs.

## Common commands

```bash
./build_lib.sh
./build_lib.sh -N -t debug -j 8
./build_lib.sh -N -t release -D CPU_ENABLE_NATIVE_TUNING=OFF
./build_lib.sh -N -D EstimationGears_for_SpaceNav_ENABLE_CUDA=ON
./build_lib.sh -N -p --gtwrap-root lib/wrap
./build_lib.sh -N -m --gtwrap-root lib/wrap
./build_lib.sh -N -i
```

Use `-D NAME=VALUE` or `-DNAME=VALUE` for repeatable CMake definitions. `-r` rebuilds an existing cache without configuring it; wrapper flags with `-r` require that the cache already enabled those wrappers.

## Cleanup safety

`--clean` accepts only conventional build paths inside this checkout. If the path already exists, its `CMakeCache.txt` must identify this exact checkout through `CMAKE_HOME_DIRECTORY`. The script will not weaken that ownership check for unusual layouts.

For disposable out-of-tree acceptance builds, invoke CMake directly and remove the directory only after separately confirming its ownership.

## Wrapper maintenance

The default is no checkout movement. These operations require explicit flags:

- `--wrap-update` updates a local wrap checkout.
- `--wrap-submodule-init` initializes the declared wrap submodule fallback.
- `--gtwrap-root PATH` selects an existing checkout without moving it.

Run `./build_lib.sh --help` for the complete current option list.
