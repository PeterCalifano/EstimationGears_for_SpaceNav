/// @file placeholder.cpp
/// @brief Implements the native scaffold entry point.

#include "placeholder.h"
#include <utils/logging/CLogger.h>

namespace placeholder
{
    void placeholder_fcn()
    {
        estimation_gears::logging::CLogger objLogger_("placeholder");
        objLogger_.setLevelFromEnvironment();
        objLogger_.info("Hello, World! I'm a placeholder function, yuppy.");
    }
} // namespace placeholder
