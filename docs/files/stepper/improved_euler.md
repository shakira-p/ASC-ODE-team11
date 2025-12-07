# Improved Euler Method (Heun's Method)

## Introduction

The Improved Euler method, also known as Heun's method, is an explicit second-order Runge-Kutta scheme. It uses a predictor-corrector approach to achieve better accuracy than the standard Euler method.

## Mathematical Overview

Given an ODE $\frac{dy}{dt} = f(t, y)$, the method computes the solution at time $t_{n+1} = t_n + \tau$ using:

$$\hat{y} = y_n + \frac{\tau}{2} f(t_n, y_n)$$

$$y_{n+1} = y_n + \tau \cdot f\left(t_n + \frac{\tau}{2}, \hat{y}\right)$$

The method first evaluates $f$ at the current point, then at a midpoint predictor $\hat{y}$, and uses this slope for the final update.

## Implementation

The `ImprovedEuler` class inherits from `TimeStepper` and performs explicit time integration.

### Constructor Parameters

- `rhs`: A `NonlinearFunction` representing $f(t, y)$

### DoStep Parameters

- `tau`: Time step size $\tau$
- `y`: Current solution vector (modified in-place to contain $y_{n+1}$)

### Key Features

- Stores intermediate values in `m_vecf` (function evaluations) and `m_y_hat` (predictor point)
- Fully explicit: no nonlinear solve required
- Second-order accurate like Crank-Nicolson, but not unconditionally stable
