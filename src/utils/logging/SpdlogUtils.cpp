#include "SpdlogUtils.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <string>

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace spdlog_utils
{
namespace
{
    constexpr auto kDefaultLoggerName = "spdlog_utils_default";
    constexpr auto kDefaultPattern = "[%H:%M:%S] [%^%l%$] [%n] %v";

    std::shared_ptr<spdlog::sinks::stderr_color_sink_mt> GetDefaultSink_()
    {
        static auto objSink_ = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
        static const bool bPatternConfigured_ = []()
        {
            objSink_->set_pattern(kDefaultPattern);
            return true;
        }();

        static_cast<void>(bPatternConfigured_);
        return objSink_;
    }

    std::shared_ptr<spdlog::logger> CreateLogger_(const std::string &charLoggerName,
                                                  const spdlog::level::level_enum enumLevel)
    {
        auto objLogger_ = std::make_shared<spdlog::logger>(charLoggerName, spdlog::sinks_init_list{GetDefaultSink_()});
        objLogger_->set_level(enumLevel);
        objLogger_->flush_on(spdlog::level::err);
        spdlog::register_logger(objLogger_);
        return objLogger_;
    }

    std::string LowerString_(std::string_view charInput)
    {
        std::string charLowered_(charInput.begin(), charInput.end());
        std::transform(charLowered_.begin(),
                       charLowered_.end(),
                       charLowered_.begin(),
                       [](unsigned char charValue_)
                       {
                           return static_cast<char>(std::tolower(charValue_));
                       });
        return charLowered_;
    }
} // namespace

void ConfigureDefaultLogging(const spdlog::level::level_enum level)
{
    auto objDefaultLogger_ = spdlog::get(kDefaultLoggerName);
    if (!objDefaultLogger_)
    {
        objDefaultLogger_ = CreateLogger_(kDefaultLoggerName, level);
    }

    objDefaultLogger_->set_level(level);
    objDefaultLogger_->flush_on(spdlog::level::err);
    spdlog::set_default_logger(objDefaultLogger_);
    spdlog::set_level(level);
}

std::shared_ptr<spdlog::logger> GetLogger(std::string_view charComponentName)
{
    if (charComponentName.empty())
    {
        ConfigureDefaultLogging();
        return spdlog::default_logger();
    }

    const std::string charLoggerName_(charComponentName);
    auto objLogger_ = spdlog::get(charLoggerName_);
    if (objLogger_)
    {
        return objLogger_;
    }

    ConfigureDefaultLogging();
    objLogger_ = CreateLogger_(charLoggerName_, spdlog::default_logger()->level());
    return objLogger_;
}

std::optional<spdlog::level::level_enum> TryParseLogLevel(std::string_view charLevelName)
{
    if (charLevelName.empty())
    {
        return std::nullopt;
    }

    int iLevel_ = 0;
    const auto *charBegin_ = charLevelName.data();
    const auto *charEnd_ = charBegin_ + charLevelName.size();
    const auto objParseResult_ = std::from_chars(charBegin_, charEnd_, iLevel_);
    if (objParseResult_.ec == std::errc{} && objParseResult_.ptr == charEnd_)
    {
        constexpr std::array<spdlog::level::level_enum, 7> arrLevelMap_{
            spdlog::level::trace,
            spdlog::level::debug,
            spdlog::level::info,
            spdlog::level::warn,
            spdlog::level::err,
            spdlog::level::critical,
            spdlog::level::off};

        if (iLevel_ >= 0 && iLevel_ < static_cast<int>(arrLevelMap_.size()))
        {
            return arrLevelMap_[static_cast<std::size_t>(iLevel_)];
        }

        return std::nullopt;
    }

    const auto charLowered_ = LowerString_(charLevelName);
    if (charLowered_ == "trace")
    {
        return spdlog::level::trace;
    }
    if (charLowered_ == "debug")
    {
        return spdlog::level::debug;
    }
    if (charLowered_ == "info")
    {
        return spdlog::level::info;
    }
    if (charLowered_ == "warn" || charLowered_ == "warning")
    {
        return spdlog::level::warn;
    }
    if (charLowered_ == "err" || charLowered_ == "error")
    {
        return spdlog::level::err;
    }
    if (charLowered_ == "critical")
    {
        return spdlog::level::critical;
    }
    if (charLowered_ == "off" || charLowered_ == "quiet")
    {
        return spdlog::level::off;
    }

    return std::nullopt;
}

bool InitializeLogLevelFromEnvironment(std::string_view charEnvVarName)
{
    const char *charEnvValue_ = std::getenv(std::string(charEnvVarName).c_str());
    if (charEnvValue_ == nullptr)
    {
        return false;
    }

    const auto objLevel_ = TryParseLogLevel(charEnvValue_);
    if (!objLevel_)
    {
        return false;
    }

    ConfigureDefaultLogging(*objLevel_);
    return true;
}
} // namespace spdlog_utils
