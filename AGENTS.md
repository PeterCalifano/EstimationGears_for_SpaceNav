# Agent instructions

## Continuity and authority

- Before context compaction, summarize the active task, Git/index state, validation evidence, and next gate in `CONTEXT.md`.
- After automatic compaction, reread `AGENTS.md` and `CONTEXT.md` before continuing.
- The current user request overrides this file. Repository architecture and verified behavior override generic donor-template defaults.
- `AGENTS.md` is normative. `CLAUDE.md` contains repository facts, commands, and architecture; do not duplicate policy there.

## Template provenance and tailoring

Repository facilities are derived from the signed `cpp_cuda_template_project` tag `v2.0.1` at commit `1d87153b2d060bf03c2c9adcd1df6c6d4f40ea09`. Treat future refreshes as a three-way import:

```text
updated project = reviewed donor implementation + preserved EstimationGears tailoring
```

Preserve these contracts:

- project/package `EstimationGears_for_SpaceNav`;
- C++ namespace `estimation_gears` and exported target `EstimationGears_for_SpaceNav::EstimationGears_for_SpaceNav`;
- C++20, Python >=3.12, MATLAB default R2023b, and CUDA >=12.6;
- optional CUDA, Python wrapping, MATLAB wrapping, TBB, OpenGL, sanitizers, documentation, install/export, and source packaging;
- no OptiX, TensorRT, ZeroMQ, profiling-script suite, or project ROS overlay;
- native builds do not recursively compose arbitrary `lib/*` checkouts;
- one Linux CPU workflow and one documentation workflow; CUDA remains a local acceptance gate;
- documentation artifacts build in CI, but the Pages `deploy` job remains disabled by `if: ${{ false }}` until the user explicitly changes release policy.

Do not import donor template-development plans, reports, cleanup scripts, issue forms, pull-request templates, ROS/CUDA workflows, or template-conformance tests.

## Git and dirty-worktree safety

- Treat all pre-existing modifications, untracked files, ignored artifacts, nested checkout changes, and staged paths as user-owned.
- Never discard, overwrite, stage, commit, amend, rebase, or push unrelated work.
- Use literal staging allowlists; never use broad `git add .` in a dirty checkout.
- Inspect the complete index with `git diff --cached` before handoff.
- Wrapper checkout updates, submodule initialization, and submodule creation are explicit maintenance operations. Ordinary configure/build commands must not move checkouts or change gitlinks.
- Do not commit or push unless the user explicitly authorizes that action.

## CMake, testing, and packaging

- `build_lib.sh` is the native library entry point. Use fresh out-of-tree CMake builds for acceptance matrices and consumers.
- `build_lib.sh --clean` may remove only a conventional in-repository build directory. An existing directory must contain `CMakeCache.txt` whose `CMAKE_HOME_DIRECTORY` resolves to this exact checkout.
- Never weaken clean-path or cache-ownership checks for unusual layouts.
- Prefer Catch2 or pytest for runtime behavior.
- Prove options, nested composition, shared/static matrices, installation, packaging, and consumers through explicit fresh commands or CI.
- Disposable consumer projects stay outside the ordinary test build.
- Add permanent CMake-script tests only when they are lightweight, target-owned, non-recursive, and cover behavior unavailable through normal targets/tests.
- Never add `VerifyTemplateProject*` or other donor self-validation tests.
- Generated Python wheels and CMake Python installs must co-locate declared non-system shared runtime targets and use loader-relative runtime paths.
- `_wrapper_build.py` is checkout-only metadata and must not be installed or included in wheels.
- CMake Python install destinations remain relative to `CMAKE_INSTALL_PREFIX`; pip owns active-environment installation.
- Configure must not write tracked source files by default.

## Python

- Use Python 3.12 or newer.
- Type hints are mandatory on every callable and meaningful variable boundary. Use precise built-in generics and avoid untyped dictionaries.
- Public functions begin with a capital letter and otherwise use snake case, for example `Load_valid_observations`.
- Classes begin with a capital letter. Class methods begin with a lowercase letter. Internal methods begin with `_`. Local variables end with `_`.
- Use Google-style module, class, method, and function docstrings.
- New classes and functions include a runnable `Example` and expected `Output`.
- Prefer dataclasses to ad-hoc dictionaries and enums to `Literal` when more than two values form a closed set.
- Prefer pathlib, context managers, explicit exceptions, deterministic resource ownership, and dependency injection over hidden global state.
- Use pytest for tests. Test public behavior and failure modes; avoid implementation-only assertions.
- Use matplotlib for general plotting and seaborn by default for statistical plots. PIL and OpenCV are appropriate for image processing.
- Use PyTorch for machine learning, supported by scikit-learn where appropriate.
- ONNX export compatibility is normally required: avoid unsupported dynamic Python behavior in model forward paths and test export plus runtime parity.

## C++ and CUDA

- C++20 is the baseline. Do not lower code to C++17 or introduce a newer requirement without an explicit compatibility decision.
- CUDA requires toolkit 12.6 or newer. Keep host/device ownership and error handling explicit.
- Follow the surrounding naming conventions and prefer classes over structs unless a type is strictly a passive aggregate.
- Prefer concepts over SFINAE and standard-library facilities over custom metaprogramming.
- Apply modern C++/Jason Turner practices:
  - use RAII and the Rule of Zero; never use naked `new`/`delete` for ownership;
  - prefer value semantics, composition, narrow scopes, and explicit lifetime boundaries;
  - make values `const` and computations `constexpr` when their semantics allow it;
  - initialize every object and avoid undefined behavior, implicit narrowing, dangling views, and unchecked ownership;
  - use strong enums/types, `std::optional` for optional values, and `[[nodiscard]]` for results callers must inspect;
  - use algorithms and ranges when they make intent clearer than index-based loops;
  - use `std::span` and `std::string_view` only when the referenced lifetime is unambiguous;
  - include what each file uses, minimize macros, and keep headers self-contained;
  - keep interfaces small, dependencies directional, and abstractions zero-overhead and measurable;
  - enable warnings and sanitizers during development; optimize only from profiler or benchmark evidence.
- Use Catch2 for unit tests. Add a failing behavioral test before new production behavior.
- Use Doxygen `@file` headers and Doxygen documentation for every public class, function, and method.
- Preserve compact grouped formatting. In multiline calls/declarations, keep the first argument beside the function name and align continuation arguments.
- Separate statements into small functional blocks. Comment purpose, invariants, and non-obvious decisions—not line-by-line syntax.
- Do not apply broad automatic formatting to staged or user-owned code.

## MATLAB

- Prefer functions unless persistent state or an object lifecycle materially improves the design. Use classes when stateful behavior is genuinely useful.
- Use `self` rather than `obj` for the instance argument.
- Function names and static class methods begin with a capital letter. Local helper functions end with `_` and are placed after the primary function; never nest function definitions.
- Use explicit, descriptive Hungarian-style datatype prefixes:
  - `d` double, `f` single, `b` logical;
  - `char` character/string data and `str` structures;
  - `ui8`/`ui16`/`ui32` and `i8`/`i16`/`i32` for integers;
  - `obj` objects, `cell` cell arrays, `table` tables, and `bus_` Simulink buses.
- Names use PascalCase after the prefix, for example `ui8MeasurementCount`. Avoid short names except in a very small local scope; temporary names should include `Tmp`.
- Use `arguments` and `arguments (Output)` blocks for public inputs and outputs where code-generation constraints allow them.
- Algorithmic MATLAB intended for deployment must remain MATLAB Coder safe. When codegen applies, keep identifiers within 31 characters and avoid dynamic constructs unsupported by Coder.
- Primary function files use this documentation structure:

```matlab
%% SIGNATURE
%
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% DD-MM-YYYY  Pietro Califano     First prototype.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
%
% -------------------------------------------------------------------------------------------------------------
```

- Keep dependencies explicit and preserve established `coder.const`, `coder.mustBeConst`, `coder.nullcopy`, and `coder.target` patterns.
- Validate numerical changes against an independent oracle where possible; same-implementation finite differences do not independently prove correctness.

## Staged-code review gate

Before presenting staged changes:

- inspect `git diff --cached --name-status`, `git diff --cached --check`, and the full `git diff --cached`;
- confirm every staged path belongs to the approved allowlist;
- confirm mixed files contain only approved hunks;
- add/update file-level and public API documentation for new or substantially modified source;
- organize related statements into readable blocks with concise intent/invariant comments;
- preserve useful existing documentation unless the change makes it false;
- search for stale donor identities, OptiX/TensorRT/ROS overlay references, generated artifacts, conflict markers, and absolute machine paths;
- rerun the applicable fresh build/test/install/wrapper/docs/package gates;
- summarize documentation/readability cleanup and any unavailable or known-red gate.

Limit cleanup to the intended facility or feature scope even when a file also contains unrelated legacy code.
