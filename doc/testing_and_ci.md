# Testing and CI

## Local native gates

```bash
cmake -S . -B /tmp/estimation-gears-test -GNinja \
  -DEstimationGears_for_SpaceNav_ENABLE_CUDA=OFF \
  -DENABLE_TESTS=ON \
  -DENABLE_PYTHON_TESTS=ON \
  -DCPU_ENABLE_NATIVE_TUNING=OFF
cmake --build /tmp/estimation-gears-test --parallel
ctest --test-dir /tmp/estimation-gears-test --output-on-failure --no-tests=error
```

Catch2 owns C++ runtime behavior, pytest owns Python behavior, and the MATLAB wrapper regression owns generated MATLAB/MEX behavior.

## Derived-project acceptance

Do not copy template-conformance verifiers into this repository. Feature matrices, nested configuration, installation, packaging, and external consumers are verified through explicit fresh commands and CI—not recursive rebuilds in the ordinary CTest suite.

Disposable consumers may prove:

- top-level and nested option isolation;
- build-tree and install-tree package use;
- shared and static builds;
- source archive hygiene;
- wrapper runtime packaging.

## Workflows

- `build_linux.yml` builds and tests the portable CPU configuration.
- `docs_pages.yml` builds Doxygen HTML/XML and uploads its artifact.
- CUDA remains a local acceptance gate rather than a GitHub workflow.
- No ROS workflow exists.

The Pages `deploy` job is intentionally present but disabled with `if: ${{ false }}`. Documentation validation remains active.
