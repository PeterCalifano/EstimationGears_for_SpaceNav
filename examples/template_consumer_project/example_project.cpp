/// @file example_project.cpp
/// @brief Demonstrates consuming an installed EstimationGears package.

#include "example_project.h"

int main()
{
    using namespace estimation_gears::logging;

    CLogger objLogger_("example_consumer_project", ELogLevel::Info);
    objLogger_.setLevelFromEnvironment();
    objLogger_.info("Hello, World! This project consumes EstimationGears through CMake.");

    // Call the placeholder function from the template_src library
    placeholder::placeholder_fcn();

    return 0;
}
