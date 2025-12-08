# Crank-Nicolson Method

## Introduction

The Crank-Nicolson method is an implicit time-stepping scheme that combines forward and backward Euler methods. It is second-order accurate and unconditionally stable, making it suitable for stiff ODEs.

## Mathematical Overview

Given an ODE $\frac{dy}{dt} = f(t, y)$, the Crank-Nicolson method computes the solution at time $t_{n+1} = t_n + \tau$ using:

$$y_{n+1} = y_n + \frac{\tau}{2} \left( f(t_n, y_n) + f(t_{n+1}, y_{n+1}) \right)$$

This is an implicit equation that requires solving the nonlinear system:

$$R(y_{n+1}) = y_{n+1} - y_n - \frac{\tau}{2} \left( f(t_n, y_n) + f(t_{n+1}, y_{n+1}) \right) = 0$$

## Implementation

The `CrankNicolson` class inherits from `TimeStepper` and uses Newton's method to solve the implicit equation.

### Constructor Parameters

- `rhs`: A `NonlinearFunction` representing $f(t, y)$

### DoStep Parameters

- `tau`: Time step size $\tau$
- `y`: Current solution vector (modified in-place to contain $y_{n+1}$)

### Key Features

- Stores $y_n$ in `m_yold` and $f(t_n, y_n)$ in `m_fold`
- Constructs the residual function $R(y_{n+1})$ using function composition
- Uses `NewtonSolver` to find $y_{n+1}$ starting from the current value
