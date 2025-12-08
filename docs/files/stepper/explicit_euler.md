# Explicit Euler

## Introduction

The explicit Euler method is the simplest and most straightforward numerical method for solving ordinary differential equations (ODEs). 
It is an explicit time-stepping scheme that updates the solution based on the current value of the derivative.

## Mathematical Overview

Given an autonomous ODE $\frac{dy}{dt} = f(y(t)) $ $\forall t \in [0,T]$, the Explicit Euler method 
computes the solution using:

$y_{i+1} = y_i + \tau f(y_i)$ $\space  0 \leq i\leq n$ with  $\tau = \frac{T}{n}$

## Implementation

The `ExplicitEuler` class inherits from `TimeStepper` and performs explicit integration.

### Constructor Parameters

- `rhs`: A `NonlinearFunction` representing $f(t, y)$

### DoStep Parameters
- `tau`: Time step size $\tau$
- `y`: Current solution vector (modified in-place to contain $y_{n+1}$)

### Key Features

The explicit euler method is the simplest of the ODE solvers implemented in the ASC-ODE package. 
It updates the solution by taking a step in the direction of the derivative at the current point.

```
this->m_rhs->evaluate(y, m_vecf);
y += tau * m_vecf;
```

## Results for Mass-Spring System

Below is a phase plot for a mass-spring system solved using the Explicit Euler method.

![Mass Spring Phase Plot Explicit Euler Method](/docs/pictures/MassSpringPhasePlot_explict.png)

![Mass Spring System Time Evolution Explicit Euler Method](/docs/pictures/MassSpringSystemTimeEvolution_explict.png)

Here we can see that the Explicit Euler method does not conserve energy, leading to a spiraling outwards trajectory in the phase plot.