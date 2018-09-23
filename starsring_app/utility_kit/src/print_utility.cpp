#include <iomanip>
#include <type_traits>

#include <utility_kit/ansi_escape_code.hpp>
#include <utility_kit/print_utility.hpp>

namespace utility {

//#############################################################################
//######################     AbstractPut             ##########################
//#############################################################################

AbstractPut::AbstractPut(PrintStyle print_style) : m_print_style(print_style) {}

std::string AbstractPut::begin_str() const {
  switch (m_print_style) {
    case PrintStyle::plain:
      return "";
    case PrintStyle::gray:
      return ColorGray;
    case PrintStyle::red:
      return ColorRed;
    case PrintStyle::green:
      return ColorGreen;
    case PrintStyle::blue:
      return ColorBlue;
    case PrintStyle::yellow:
      return ColorYellow;
    default:
      return "";
  }
}

std::string AbstractPut::end_str() const {
  if (m_print_style != PrintStyle::plain) return ColorEnd;
  return "";
}

//#############################################################################
//######################     ValuePut                ##########################
//#############################################################################

template <typename T>
ValuePut<T>::ValuePut(std::string title, unsigned title_width, T value,
                      unsigned value_width, unsigned value_precision,
                      std::string unit)
    : AbstractPut(PrintStyle::plain),
      m_title(title),
      m_title_width(title_width),
      m_value(value),
      m_value_width(value_width),
      m_value_precision(value_precision),
      m_unit(unit) {}

template <typename T>
ValuePut<T>::ValuePut(PrintStyle print_style, std::string title,
                      unsigned title_width, T value, unsigned value_width,
                      unsigned value_precision, std::string unit)
    : AbstractPut(print_style),
      m_title(title),
      m_title_width(title_width),
      m_value(value),
      m_value_width(value_width),
      m_value_precision(value_precision),
      m_unit(unit) {}

template <typename T>
std::ostream &operator<<(std::ostream &stream, const ValuePut<T> &value_put) {
  std::ios::fmtflags f(stream.flags());
  stream << value_put.begin_str();
  if (std::is_arithmetic<T>::value)
    stream << std::fixed << std::boolalpha
           << std::setprecision(value_put.m_value_precision) << std::showpos;
  stream << std::left << std::setw(value_put.m_title_width) << value_put.m_title
         << " : ";
  stream << (std::is_arithmetic<T>::value ? std::right : std::left);
  stream << std::setw(value_put.m_value_width) << value_put.m_value << " "
         << value_put.m_unit;
  stream << value_put.end_str();
  stream << std::endl;
  stream.flags(f);
  return stream;
}  // end of std::ostream & operator<<(std::ostream & stream, const ValuePut<T>
   // & value_put)

template class ValuePut<bool>;
template class ValuePut<int>;
template class ValuePut<unsigned>;
template class ValuePut<double>;
template class ValuePut<std::string>;
template std::ostream &operator<<(std::ostream &, const ValuePut<bool> &);
template std::ostream &operator<<(std::ostream &, const ValuePut<int> &);
template std::ostream &operator<<(std::ostream &, const ValuePut<unsigned> &);
template std::ostream &operator<<(std::ostream &, const ValuePut<double> &);
template std::ostream &operator<<(std::ostream &,
                                  const ValuePut<std::string> &);

//#############################################################################
//######################    TernaryPut               ##########################
//#############################################################################

TernaryPut::TernaryPut(std::string title, unsigned title_width, double x,
                       double y, double z, unsigned value_width,
                       unsigned value_precision, std::string unit)
    : AbstractPut(PrintStyle::plain),
      m_title(title),
      m_title_width(title_width),
      m_x(x),
      m_y(y),
      m_z(z),
      m_value_width(value_width),
      m_value_precision(value_precision),
      m_unit(unit) {}

TernaryPut::TernaryPut(PrintStyle print_style, std::string title,
                       unsigned title_width, double x, double y, double z,
                       unsigned value_width, unsigned value_precision,
                       std::string unit)
    : AbstractPut(print_style),
      m_title(title),
      m_title_width(title_width),
      m_x(x),
      m_y(y),
      m_z(z),
      m_value_width(value_width),
      m_value_precision(value_precision),
      m_unit(unit) {}

TernaryPut::TernaryPut(std::string title, unsigned title_width, arma::vec3 v,
                       unsigned value_width, unsigned value_precision,
                       std::string unit)
    : AbstractPut(PrintStyle::plain),
      m_title(title),
      m_title_width(title_width),
      m_x(v(0)),
      m_y(v(1)),
      m_z(v(2)),
      m_value_width(value_width),
      m_value_precision(value_precision),
      m_unit(unit) {}

TernaryPut::TernaryPut(PrintStyle print_style, std::string title,
                       unsigned title_width, arma::vec3 v, unsigned value_width,
                       unsigned value_precision, std::string unit)
    : AbstractPut(print_style),
      m_title(title),
      m_title_width(title_width),
      m_x(v(0)),
      m_y(v(1)),
      m_z(v(2)),
      m_value_width(value_width),
      m_value_precision(value_precision),
      m_unit(unit) {}

std::ostream &operator<<(std::ostream &stream, const TernaryPut &ternary_put) {
  std::ios::fmtflags f(stream.flags());
  stream << ternary_put.begin_str();
  stream << std::fixed << std::setprecision(ternary_put.m_value_precision)
         << std::showpos;
  stream << std::left << std::setw(ternary_put.m_title_width)
         << ternary_put.m_title << " : " << std::right
         << std::setw(ternary_put.m_value_width) << ternary_put.m_x << " "
         << std::setw(ternary_put.m_value_width) << ternary_put.m_y << " "
         << std::setw(ternary_put.m_value_width) << ternary_put.m_z << " "
         << ternary_put.m_unit;
  stream << ternary_put.end_str();
  stream << std::endl;
  stream.flags(f);
  return stream;
}  // end of std::ostream & operator<<(std::ostream & stream, const TernaryPut &
   // ternary_put)

}  // end of namespace utility
