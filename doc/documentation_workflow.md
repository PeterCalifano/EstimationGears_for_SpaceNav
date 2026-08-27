# Documentation Workflow

Doxygen reads `README.md`, `src/`, and `doc/`. It excludes submodules, build trees, development plans, and historical reports.

## Local build

```bash
cmake --preset docs
cmake --build --preset docs
```

Generated HTML is written under `build_docs/doc/html/` and XML under `build_docs/doc/xml/`.

The direct equivalent is:

```bash
cmake -S . -B build_docs -GNinja \
  -DEstimationGears_for_SpaceNav_ENABLE_CUDA=OFF \
  -DENABLE_TESTS=OFF \
  -DBUILD_DOC_HTML=ON \
  -DBUILD_DOC_XML=ON
cmake --build build_docs --target doc
```

## GitHub Pages workflow

`.github/workflows/docs_pages.yml` builds documentation and uploads the Pages artifact on relevant pushes, pull requests, and manual dispatches. The `deploy` job and its permissions remain defined for later activation, but deployment is currently disabled by the exact job condition:

```yaml
if: ${{ false }}
```

Changing that condition is a separate release-policy decision.
