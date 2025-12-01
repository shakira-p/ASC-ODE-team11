#include <iostream>
#include <cmath>
#include "nonlinfunc.hpp"
#include "autodiff.hpp"
#include "vector.hpp"
#include "matrix.hpp"

using namespace ASC_ode;
using namespace nanoblas;

int main()
{
  std::cout << "Testing PendulumAD with Automatic Differentiation\n";
  std::cout << "=================================================\n\n";

  double length = 1.0;      // 1 meter
  double gravity = 9.81;    // m/s^2

  PendulumAD pendulum(length, gravity);

  // Test point: alpha = pi/4, alpha' = 0.5 rad/s
  Vector<double> x(2);
  x(0) = M_PI / 4.0;  // angular displacement
  x(1) = 0.5;         // angular velocity

  std::cout << "Pendulum parameters:\n";
  std::cout << "  Length l = " << length << " m\n";
  std::cout << "  Gravity g = " << gravity << " m/s^2\n\n";

  std::cout << "Test point:\n";
  std::cout << "  alpha = " << x(0) << " rad (" << x(0) * 180.0 / M_PI << " degrees)\n";
  std::cout << "  alpha' = " << x(1) << " rad/s\n\n";

  Vector<double> f(2);
  pendulum.evaluate(x, f);

  std::cout << "Function evaluation f(x):\n";
  std::cout << "  f[0] = alpha' = " << f(0) << "\n";
  std::cout << "  f[1] = -(g/l)*sin(alpha) = " << f(1) << "\n\n";

  // Evaluate Jacobian using AutoDiff
  Matrix<double> jacobian(2, 2);
  pendulum.evaluateDeriv(x, jacobian);

  std::cout << "Jacobian df/dx computed using AutoDiff:\n";
  std::cout << "  df[0]/dx[0] = " << jacobian(0, 0) << "\n";
  std::cout << "  df[0]/dx[1] = " << jacobian(0, 1) << "\n";
  std::cout << "  df[1]/dx[0] = " << jacobian(1, 0) << "\n";
  std::cout << "  df[1]/dx[1] = " << jacobian(1, 1) << "\n\n";

  std::cout << "Expected Jacobian (analytical):\n";
  std::cout << "  df[0]/dx[0] = 0\n";
  std::cout << "  df[0]/dx[1] = 1\n";
  std::cout << "  df[1]/dx[0] = -(g/l)*cos(alpha) = " << -(gravity/length) * std::cos(x(0)) << "\n";
  std::cout << "  df[1]/dx[1] = 0\n\n";

  return 0;
}
