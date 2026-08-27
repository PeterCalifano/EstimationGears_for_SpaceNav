# Dependency-free logging

`estimation_gears::logging::CLogger` is a small C++20 logger owned by the project. It has no third-party logging dependency.

## Contract

- Levels are `Quiet`, `Critical`, `Error`, `Warning`, `Info`, `Debug`, and `Trace`.
- Critical, error, and warning messages use the diagnostic stream.
- Info, debug, and trace messages use the ordinary output stream.
- Each complete line is assembled before a process-wide output lock is acquired.
- ANSI color is explicit and disabled by default.
- `ESTIMATION_GEARS_LOG_LEVEL` configures the default environment-driven level.

## Example

```cpp
#include <utils/logging/CLogger.h>

int main()
{
    using namespace estimation_gears::logging;

    CLogger objLogger_("navigation", ELogLevel::Info);
    objLogger_.setLevelFromEnvironment();
    objLogger_.info("filter initialized");
}
```

Output:

```text
[navigation][INFO] filter initialized
```

The logger is non-copyable and non-movable because it stores stream references. Callers supplying custom streams must keep them alive for the logger lifetime.
