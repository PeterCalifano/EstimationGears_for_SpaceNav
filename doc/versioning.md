# Versioning Guide

## Resolution order

Both CMake and `generate_version.sh` resolve version information in this order:

1. the nearest semantic Git tag and `git describe` metadata;
2. an existing `VERSION` file;
3. the hardcoded `0.3.0` fallback.

The build-tree `VERSION` file is always generated. Source-tree generation is opt-in:

```bash
cmake -S . -B /tmp/estimation-gears-version \
  -DWRITE_SOURCE_VERSION_FILE=ON
```

Ordinary configure leaves the checkout untouched.

## Version representation

`PROJECT_VERSION` contains the semantic core accepted by CMake. `FULL_VERSION` may append prerelease and Git metadata, for example:

```text
0.3.0+91.gfe517cb.dirty
```

The configured `config.h` exposes the core and full versions. Python package metadata uses the normalized package version produced by the wrapper facility.

## Release artifacts

Release tags use `vX.Y.Z`. CPack creates binary and source TGZ archives named from the resolved full version. Source packaging excludes Git state, Codex context, builds, generated MATLAB/C++ artifacts, Python metadata, native wrapper binaries, and `_wrapper_build.py`.

This repository has no ROS package metadata synchronization.
