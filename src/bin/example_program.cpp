#include <template_src/placeholder.h>

#ifndef SPDLOG_UTILS_ENABLED
#define SPDLOG_UTILS_ENABLED 0
#endif

#if SPDLOG_UTILS_ENABLED
#include <utils/logging/SpdlogUtils.h>
#endif

int main()
{
#if SPDLOG_UTILS_ENABLED
    if (!spdlog_utils::InitializeLogLevelFromEnvironment())
    {
        spdlog_utils::ConfigureDefaultLogging();
    }

    auto objLogger_ = spdlog_utils::GetLogger("example_program");
    objLogger_->info("Running template example program.");
#endif

    placeholder::placeholder_fcn();
    return 0;
}
