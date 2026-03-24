#include <catch2/catch_test_macros.hpp>
#include <utils/logging/SpdlogUtils.h>

#include <cstdlib>
#include <optional>
#include <string>

#include <spdlog/spdlog.h>

namespace
{
class CEnvironmentVariableGuard
{
  public:
    explicit CEnvironmentVariableGuard(const std::string &charVariableName)
        : charVariableName_(charVariableName)
    {
        const char *charExistingValue_ = std::getenv(charVariableName_.c_str());
        if (charExistingValue_ != nullptr)
        {
            charPreviousValue_ = std::string(charExistingValue_);
        }
    }

    ~CEnvironmentVariableGuard()
    {
        if (charPreviousValue_.has_value())
        {
            setenv(charVariableName_.c_str(), charPreviousValue_->c_str(), 1);
        }
        else
        {
            unsetenv(charVariableName_.c_str());
        }
    }

    void setValue(const std::string &charValue) const
    {
        setenv(charVariableName_.c_str(), charValue.c_str(), 1);
    }

  private:
    std::string charVariableName_{};
    std::optional<std::string> charPreviousValue_{};
};
} // namespace

TEST_CASE("logging level parsing handles supported inputs", "[logging]")
{
    using spdlog_utils::TryParseLogLevel;

    REQUIRE(TryParseLogLevel("trace") == spdlog::level::trace);
    REQUIRE(TryParseLogLevel("DEBUG") == spdlog::level::debug);
    REQUIRE(TryParseLogLevel("warning") == spdlog::level::warn);
    REQUIRE(TryParseLogLevel("error") == spdlog::level::err);
    REQUIRE(TryParseLogLevel("6") == spdlog::level::off);
    REQUIRE_FALSE(TryParseLogLevel("verbose").has_value());
}

TEST_CASE("environment initialization applies parsed log level", "[logging]")
{
    CEnvironmentVariableGuard objEnvGuard_("SPDLOG_LEVEL");
    objEnvGuard_.setValue("debug");

    REQUIRE(spdlog_utils::InitializeLogLevelFromEnvironment());
    REQUIRE(spdlog::default_logger() != nullptr);
    REQUIRE(spdlog::default_logger()->level() == spdlog::level::debug);
}

TEST_CASE("logger lookup reuses registered logger instances", "[logging]")
{
    spdlog_utils::ConfigureDefaultLogging(spdlog::level::info);

    auto objLoggerA_ = spdlog_utils::GetLogger("logging_test_component");
    auto objLoggerB_ = spdlog_utils::GetLogger("logging_test_component");

    REQUIRE(objLoggerA_ != nullptr);
    REQUIRE(objLoggerA_.get() == objLoggerB_.get());
}

TEST_CASE("default logging configuration is safe to repeat", "[logging]")
{
    spdlog_utils::ConfigureDefaultLogging(spdlog::level::warn);
    auto objFirstDefaultLogger_ = spdlog::default_logger();

    spdlog_utils::ConfigureDefaultLogging(spdlog::level::debug);
    auto objSecondDefaultLogger_ = spdlog::default_logger();

    REQUIRE(objFirstDefaultLogger_ != nullptr);
    REQUIRE(objFirstDefaultLogger_.get() == objSecondDefaultLogger_.get());
    REQUIRE(objSecondDefaultLogger_->level() == spdlog::level::debug);
}
