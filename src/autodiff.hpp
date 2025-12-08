#ifndef AUTODIFF_HPP
#define AUTODIFF_HPP

#include <cstddef>
#include <ostream>
#include <cmath>
#include <array>
#include <vector>
#include <algorithm>
#include <type_traits>


namespace ASC_ode
{

  template <typename T = double>
  auto derivative (T v, size_t /*index*/) { return T(0); } 


  // Dynamic AutoDiff class - size determined at runtime
  template <typename T = double>
  class AutoDiff
  {
  private:
    T m_val;
    std::vector<T> m_deriv;
  public:
    // Default constructor
    AutoDiff () : m_val(0), m_deriv() {}

    // Construct from value with specified size (all derivatives zero)
    AutoDiff (T v, size_t n = 0) : m_val(v), m_deriv(n, T(0)) {}

    // Construct as variable: value v, derivative 1 at index derivIndex
    AutoDiff(T val, size_t derivIndex, size_t n) : m_val(val), m_deriv(n, T(0))
    {
      if (derivIndex < n)
        m_deriv[derivIndex] = T(1);
    }

    T value() const { return m_val; }
    std::vector<T>& deriv() { return m_deriv; }
    const std::vector<T>& deriv() const { return m_deriv; }
    size_t size() const { return m_deriv.size(); }

    // Resize the derivative vector
    void resize(size_t n) { m_deriv.resize(n, T(0)); }
  };


  template <typename T = double>
  auto derivative (const AutoDiff<T>& v, size_t index)
  {
    if (index < v.deriv().size())
      return v.deriv()[index];
    return T(0);
  }


  template <typename T>
  std::ostream & operator<< (std::ostream& os, const AutoDiff<T>& ad)
  {
    os << "Value: " << ad.value() << ", Deriv: [";
    for (size_t i = 0; i < ad.size(); i++)
    {
      os << ad.deriv()[i];
      if (i < ad.size() - 1) os << ", ";
    }
    os << "]";
    return os;
  }

  // Helper to get max size of two AutoDiff objects
  template <typename T>
  size_t maxSize(const AutoDiff<T>& a, const AutoDiff<T>& b)
  {
    return std::max(a.size(), b.size());
  }

  // Helper to safely get derivative (returns 0 if index out of range)
  template <typename T>
  T safeGetDeriv(const AutoDiff<T>& a, size_t i)
  {
    return (i < a.size()) ? a.deriv()[i] : T(0);
  }

  template <typename T = double>
  AutoDiff<T> operator+ (const AutoDiff<T>& a, const AutoDiff<T>& b)
  {
    size_t n = maxSize(a, b);
    AutoDiff<T> result(a.value() + b.value(), n);
    for (size_t i = 0; i < n; i++)
      result.deriv()[i] = safeGetDeriv(a, i) + safeGetDeriv(b, i);
    return result;
  }

  template <typename T = double>
  auto operator+ (T a, const AutoDiff<T>& b) { return AutoDiff<T>(a, b.size()) + b; }

  template <typename T = double>
  auto operator+ (const AutoDiff<T>& a, T b) { return a + AutoDiff<T>(b, a.size()); }

  template <typename T = double>
  AutoDiff<T> operator* (const AutoDiff<T>& a, const AutoDiff<T>& b)
  {
    size_t n = maxSize(a, b);
    AutoDiff<T> result(a.value() * b.value(), n);
    for (size_t i = 0; i < n; i++)
      result.deriv()[i] = safeGetDeriv(a, i) * b.value() + a.value() * safeGetDeriv(b, i);
    return result;
  }

  template <typename T = double>
  auto operator* (T s, const AutoDiff<T>& b) { return AutoDiff<T>(s, b.size()) * b; }

  template <typename T = double>
  auto operator* (const AutoDiff<T>& a, T s) { return a * AutoDiff<T>(s, a.size()); }

  // unary minus: -ad
  template <typename T = double>
  AutoDiff<T> operator- (const AutoDiff<T>& a)
  {
    AutoDiff<T> result(-a.value(), a.size());
    for (size_t i = 0; i < a.size(); ++i)
      result.deriv()[i] = -a.deriv()[i];
    return result;
  }

  // binary subtraction: ad1 - ad2
  template <typename T = double>
  AutoDiff<T> operator- (const AutoDiff<T>& a, const AutoDiff<T>& b)
  {
    size_t n = maxSize(a, b);
    AutoDiff<T> result(a.value() - b.value(), n);
    for (size_t i = 0; i < n; ++i)
      result.deriv()[i] = safeGetDeriv(a, i) - safeGetDeriv(b, i);
    return result;
  }

  // AutoDiff - scalar
  template <typename T = double>
  AutoDiff<T> operator- (const AutoDiff<T>& a, T b)
  {
    return a - AutoDiff<T>(b, a.size());
  }

  // scalar - AutoDiff
  template <typename T = double>
  AutoDiff<T> operator- (T a, const AutoDiff<T>& b)
  {
    return AutoDiff<T>(a, b.size()) - b;
  }

  // division: ad1 / ad2
  template <typename T = double>
  AutoDiff<T> operator/ (const AutoDiff<T>& a, const AutoDiff<T>& b)
  {
    size_t n = maxSize(a, b);
    AutoDiff<T> result(a.value() / b.value(), n);
    T denom = b.value() * b.value();
    for (size_t i = 0; i < n; ++i)
      result.deriv()[i] = (safeGetDeriv(a, i) * b.value() - a.value() * safeGetDeriv(b, i)) / denom;
    return result;
  }

  // AutoDiff / scalar
  template <typename T = double>
  AutoDiff<T> operator/ (const AutoDiff<T>& a, T b)
  {
    return a / AutoDiff<T>(b, a.size());
  }

  // scalar / AutoDiff
  template <typename T = double>
  auto operator/ (T a, const AutoDiff<T>& b) { return AutoDiff<T>(a, b.size()) / b; }

  // operator +=
  template <typename T = double>
  auto operator+= (AutoDiff<T>& a, const AutoDiff<T>& b)
  {
    a = a + b;
    return a;
  }

  template <typename T = double>
  auto operator-= (AutoDiff<T>& a, const AutoDiff<T>& b)
  {
    a = a - b;
    return a;
  }

  template <typename T = double>
  auto operator*= (AutoDiff<T>& a, const AutoDiff<T>& b)
  {
    a = a * b;
    return a;
  }

  template <typename T = double>
  auto operator/= (AutoDiff<T>& a, const AutoDiff<T>& b)
  {
    a = a / b;
    return a;
  }

  using std::sin;
  using std::cos;

  template <typename T = double>
  AutoDiff<T> sin(const AutoDiff<T> &a)
  {
    AutoDiff<T> result(sin(a.value()), a.size());
    for (size_t i = 0; i < a.size(); i++)
      result.deriv()[i] = cos(a.value()) * a.deriv()[i];
    return result;
  }

  template <typename T = double>
  AutoDiff<T> cos(const AutoDiff<T> &a)
  {
    AutoDiff<T> result(cos(a.value()), a.size());
    for (size_t i = 0; i < a.size(); i++)
      result.deriv()[i] = -sin(a.value()) * a.deriv()[i];
    return result;
  }

  using std::exp;
  using std::log;

  template <typename T = double>
  AutoDiff<T> exp(const AutoDiff<T>& a)
  {
    T v = exp(a.value());
    AutoDiff<T> result(v, a.size());
    for (size_t i = 0; i < a.size(); ++i)
      result.deriv()[i] = v * a.deriv()[i]; // d/dx exp(u) = exp(u)*u'
    return result;
  }

  template <typename T = double>
  AutoDiff<T> log(const AutoDiff<T>& a)
  {
    AutoDiff<T> result(log(a.value()), a.size());
    for (size_t i = 0; i < a.size(); ++i)
      result.deriv()[i] = a.deriv()[i] / a.value(); // d/dx log(u) = u'/u
    return result;
  }

  using std::sqrt;

  template <typename T = double>
  AutoDiff<T> sqrt(const AutoDiff<T>& a)
  {
    T v = sqrt(a.value());
    AutoDiff<T> result(v, a.size());
    for (size_t i = 0; i < a.size(); ++i)
      result.deriv()[i] = a.deriv()[i] / (T(2) * v);  // d/dx sqrt(u) = u'/(2*sqrt(u))
    return result;
  }

  // norm2 for scalar AutoDiff (squared magnitude, used by nanoblas internally)
  template <typename T = double>
  AutoDiff<T> norm2(const AutoDiff<T>& a)
  {
    return a * a;  // For real numbers, |x|^2 = x^2
  }

  // Vector norm for any indexable container of AutoDiff (or double)
  // Works with nanoblas Vec, std::array, std::vector, etc.
  template <typename VEC>
  auto vecNorm(const VEC& v) -> std::remove_const_t<std::remove_reference_t<decltype(v(0))>>
  {
    using T = std::remove_const_t<std::remove_reference_t<decltype(v(0))>>;
    T sum = T(0.0);
    for (size_t i = 0; i < v.size(); i++)
      sum = sum + v(i) * v(i);
    return sqrt(sum);
  }


} // namespace ASC_ode

#endif
