#ifndef UTILITY_KIT_PRINT_UTILITY_HPP
#define UTILITY_KIT_PRINT_UTILITY_HPP

#include <armadillo>
#include <iostream>
#include <string>

#include "utility_kit/ansi_escape_code.hpp"

//#############################################################################
//######################     PREPROCESOR PRINT       ##########################
//#############################################################################

#define PRINT_ANGLE(phi)                                              \
  std::cout << #phi << " = " << (phi)                                 \
            << " rad =  " << (phi) * (180 / arma::datum::pi) << " st" \
            << std::endl;
#define PRINT(x) std::cout << #x << " = " << (x) << std::endl;

namespace utility {

//#############################################################################
//######################     KOLORKI                 ##########################
//#############################################################################

enum class PrintStyle {
  plain,
  gray,
  red,
  green,
  blue,
  yellow,
};

//#############################################################################
//######################     AbstractPut             ##########################
//#############################################################################

class AbstractPut {
 public:
  AbstractPut(PrintStyle print_style);
  const PrintStyle m_print_style;
  std::string begin_str() const;
  std::string end_str() const;
  virtual ~AbstractPut() = default;
};

//#############################################################################
//######################     ValuePut                ##########################
//#############################################################################

template <typename T>
class ValuePut : private AbstractPut {
 private:
  const std::string m_title;
  const unsigned m_title_width;
  const T m_value;
  const unsigned m_value_width;
  const unsigned m_value_precision;
  const std::string m_unit;

 public:
  ValuePut(std::string title, unsigned title_width, T value,
           unsigned value_width = 0, unsigned value_precision = 6,
           std::string unit = std::string());
  ValuePut(PrintStyle print_style, std::string title, unsigned title_width,
           T value, unsigned value_width = 0, unsigned value_precision = 6,
           std::string unit = std::string());
  template <typename U>
  friend std::ostream& operator<<(std::ostream& stream,
                                  const ValuePut<U>& ternary_put);
};

//#############################################################################
//######################    TernaryPut               ##########################
//#############################################################################

class TernaryPut : private AbstractPut {
 private:
  const std::string m_title;
  const unsigned m_title_width;
  const double m_x;
  const double m_y;
  const double m_z;
  const unsigned m_value_width;
  const unsigned m_value_precision;
  const std::string m_unit;

 public:
  TernaryPut(std::string title, unsigned title_width, double x, double y,
             double z, unsigned value_width = 0, unsigned value_precision = 6,
             std::string unit = std::string());
  TernaryPut(PrintStyle print_style, std::string title, unsigned title_width,
             double x, double y, double z, unsigned value_width = 0,
             unsigned value_precision = 6, std::string unit = std::string());
  TernaryPut(std::string title, unsigned title_width, arma::vec3 v,
             unsigned value_width = 0, unsigned value_precision = 6,
             std::string unit = std::string());
  TernaryPut(PrintStyle print_style, std::string title, unsigned title_width,
             arma::vec3 v, unsigned value_width = 0,
             unsigned value_precision = 6, std::string unit = std::string());
  friend std::ostream& operator<<(std::ostream& stream,
                                  const TernaryPut& ternary_put);
};

}  // end of namespace utility

#endif  // UTILITY_KIT_PRINT_UTILITY_HPP
