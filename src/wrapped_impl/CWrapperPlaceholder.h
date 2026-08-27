/// @file CWrapperPlaceholder.h
/// @brief Declares the minimal API used to validate generated wrappers.

#pragma once

#include <cstdint>
#include <string>

namespace estimation_gears
{
    /// @brief Small value-owning class exercised by Python and MATLAB wrappers.
    class CWrapperPlaceholder
    {
      public:
        /// @brief Construct the placeholder with deterministic default values.
        CWrapperPlaceholder() = default;

        /// @brief Return the stored numeric value.
        double getDataMember() const;
        /// @brief Replace the stored numeric value.
        /// @param value New numeric value.
        void setDataMember(double value);
        /// @brief Return the stored text value.
        std::string getTextData() const;
        /// @brief Replace the stored text from a constant reference.
        /// @param charValue New text value.
        void setTextDataByConstRef(const std::string &charValue);
        /// @brief Replace the stored text by value, permitting move construction.
        /// @param charValue New text value.
        void setTextDataByValue(std::string charValue);
        /// @brief Round-trip an unsigned frame identifier through the wrapper.
        /// @param ui32FrameId Frame identifier.
        /// @return The unchanged identifier.
        std::uint32_t echoFrameId(std::uint32_t ui32FrameId) const;
        /// @brief Print the current text value to standard output.
        void printToStdout() const;

        /// @brief Multiply a value by two.
        /// @param value Input value.
        /// @return Twice the input value.
        static double multiplyBy2(double value);

      private:
        double a_float_number_{0.0};
        std::string charTextData_{"initial"};
    };

} // namespace estimation_gears
