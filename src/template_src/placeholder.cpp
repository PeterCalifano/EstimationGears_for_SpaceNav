#include "placeholder.h"

#ifndef SPDLOG_UTILS_ENABLED
#define SPDLOG_UTILS_ENABLED 0
#endif

#if SPDLOG_UTILS_ENABLED
#include <utils/logging/SpdlogUtils.h>
#endif

namespace placeholder 
{
    void placeholder_fcn()
    {
#if SPDLOG_UTILS_ENABLED
        auto objLogger_ = spdlog_utils::GetLogger("placeholder");
        objLogger_->info("Hello, World! I'm a placeholder function, yuppy.");
#else
        std::cout << "Hello, World! I'm a placeholder function, yuppy." << std::endl;
#endif
    }
}
