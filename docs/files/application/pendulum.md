# Pendulum Automatic Differentiation

## Introduction

The `PendulumAD` class models the differential equations of a simple pendulum using automatic differentiation (AD) to compute derivatives. It implements the `NonlinearFunction` interface for use with time-stepping and Newton solvers.

## Mathematical Overview

A simple pendulum with length $L$ under gravity $g$ is governed by:

$$\frac{d\theta}{dt} = \omega$$

$$\frac{d\omega}{dt} = -\frac{g}{L}\sin(\theta)$$

where $\theta$ is the angular displacement and $\omega$ is the angular velocity. The state vector is $x = [\theta, \omega]^T$ and the right-hand side is:

$$F(x) = \begin{bmatrix} \omega \\ -\frac{g}{L}\sin(\theta) \end{bmatrix}$$

The Jacobian matrix is computed automatically using automatic differentiation:

$$J_F(x) = \begin{bmatrix} 0 & 1 \\ -\frac{g}{L}\cos(\theta) & 0 \end{bmatrix}$$

## Implementation

The `PendulumAD` class uses template-based automatic differentiation to compute both function values and derivatives.

### Constructor Parameters

- `length`: Length $L$ of the pendulum (in meters)
- `gravity`: Gravitational acceleration $g$ (default: 9.81 m/s²)

### Key Features

- **AD**: Stands for **Automatic Differentiation**, a technique that computes exact derivatives using dual numbers
- Uses the `AutoDiff<2>` type to track derivatives with respect to both state variables
- The `T_evaluate` template method works with both `double` (for function evaluation) and `AutoDiff<2>` (for derivative computation)
- `evaluate()`: Computes $F(x)$ using standard doubles
- `evaluateDeriv()`: Computes the Jacobian $J_F(x)$ by evaluating with `AutoDiff<2>` types and extracting derivatives
