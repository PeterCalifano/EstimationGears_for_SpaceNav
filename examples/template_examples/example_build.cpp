/// @file example_build.cpp
/// @brief Demonstrates EstimationGears logging and placeholder usage.

#include <template_src/placeholder.h>
#include <utils/logging/CLogger.h>

int main()
{
    using namespace estimation_gears::logging;

    CLogger objLogger_("example_build", ELogLevel::Info);
    objLogger_.setLevelFromEnvironment();
    objLogger_.info("Hello, World! This is an EstimationGears example.");
    placeholder::placeholder_fcn();

    // Example output:
    // [example_build][INFO] Hello, World! This is an EstimationGears example.

    return 0;
}
