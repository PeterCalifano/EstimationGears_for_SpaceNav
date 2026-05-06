#pragma once

namespace estimation_gears
{
class CWrapperPlaceholder
{
  public:
    CWrapperPlaceholder() = default;

    double getDataMember() const;
    void setDataMember(double value);

    static double multiplyBy2(double value);

  private:
    double a_float_number_{0.0};
};

} // namespace estimation_gears
