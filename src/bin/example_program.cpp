/// @file example_program.cpp
/// @brief Runs the EstimationGears native scaffold example.

#include <template_src/placeholder.h>
#include <utils/logging/CLogger.h>

int main()
{
    estimation_gears::logging::CLogger objLogger_("example_program");
    objLogger_.setLevelFromEnvironment();
    objLogger_.info("Running EstimationGears example program.");
    objLogger_.debug("Detailed diagnostics are enabled.");

    placeholder::placeholder_fcn();
    return 0;
}
