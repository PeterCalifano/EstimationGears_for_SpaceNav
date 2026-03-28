#pragma once

#include <memory>
#include <optional>
#include <string_view>

#include <spdlog/common.h>
#include <spdlog/logger.h>

namespace spdlog_utils
{
    /**
     * @brief Configures default logging with spdlog library
     * 
     * @param level The logging level to set (default: info)
     */
    void ConfigureDefaultLogging(spdlog::level::level_enum level = spdlog::level::info);

    /**
     * @brief Get the Logger object for a specific component
     * 
     * @param charComponentName The name of the component to get logger for
     * @return std::shared_ptr<spdlog::logger> Shared pointer to the logger instance
     */
    std::shared_ptr<spdlog::logger> GetLogger(std::string_view charComponentName);

    /**
     * @brief Attempts to parse a log level string to spdlog::level::level_enum
     * 
     * @param charLevelName The log level name as a string (e.g., "debug", "info", "warn")
     * @return std::optional<spdlog::level::level_enum> The parsed log level if successful, std::nullopt otherwise
     */
    std::optional<spdlog::level::level_enum> TryParseLogLevel(std::string_view charLevelName);

    /**
     * @brief Initializes the log level from an environment variable
     * 
     * @param charEnvVarName The environment variable name to read (default: "SPDLOG_LEVEL")
     * @return true if log level was successfully initialized from environment variable
     * @return false if environment variable was not set or log level parsing failed
     */
    bool InitializeLogLevelFromEnvironment(std::string_view charEnvVarName = "SPDLOG_LEVEL");

} // namespace spdlog_utils
