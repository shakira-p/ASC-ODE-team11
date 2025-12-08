# Automatic Differentiation (AutoDiff)

## Overview

The `AutoDiff<T>` class implements **forward-mode automatic differentiation**, computing exact derivatives alongside function values in a single pass.

## Why Automatic Differentiation?

- **Symbolic Differentiation**: Exact results, but complex expressions and slow evaluation
- **Finite Differences**: Simple to implement, but suffers from numerical errors and requires O(n) evaluations
- **AutoDiff**: Exact derivatives with efficient computation, requires operator overloading

AutoDiff propagates derivatives through the chain rule automatically, avoiding truncation errors inherent in finite differences.

## Dynamic Vector Sizing

The derivative vector uses `std::vector<T>` rather than fixed-size arrays, enabling:
- Runtime-determined variable counts (essential for Python bindings via pybind11)
- Flexible problem sizes without recompilation
- Seamless interoperability with dynamic Python arrays

## Usage

```cpp
using namespace ASC_ode;

// Create variable x = 2.0, derivative index 0, total 2 variables
AutoDiff<double> x(2.0, 0, 2);
// Create variable y = 3.0, derivative index 1, total 2 variables
AutoDiff<double> y(3.0, 1, 2);

auto f = x * y + sin(x);  // f = x*y + sin(x)
// f.value() = 6.0 + sin(2.0)
// f.deriv()[0] = y + cos(x) = 3 + cos(2)  (df/dx)
// f.deriv()[1] = x = 2                     (df/dy)
```

## Supported Operators

- `+`, `-`, `*`, `/` - Binary arithmetic (AutoDiff <-> AutoDiff, AutoDiff <-> scalar)
- `-` (unary) - Negation
- `+=`, `-=`, `*=`, `/=` - Compound assignment

## Supported Functions

- `sin(u)` - Derivative: cos(u) * u'
- `cos(u)` - Derivative: -sin(u) * u'
- `exp(u)` - Derivative: exp(u) * u'
- `log(u)` - Derivative: u' / u
- `sqrt(u)` - Derivative: u' / (2 * sqrt(u))
- `norm2(u)` - Returns u^2 (squared magnitude)
- `vecNorm(v)` - Returns sqrt(sum(vi^2)) for vectors

## Helper Functions

- `derivative(ad, index)` - Extract derivative at specific index
- `maxSize(a, b)` - Get maximum derivative vector size
- `safeGetDeriv(a, i)` - Safe access returning 0 if out of range

