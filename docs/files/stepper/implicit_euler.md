# Implicit Euler

## Introduction

The implicit Euler method, also known as the backward Euler method,
is an implicit first-order time-stepping scheme for solving ordinary differential equations (ODEs).
It is unconditionally stable and is particularly useful for stiff ODEs.

## Mathematical Overview

For every timestep for the value $y^{n+1}$ we need to solve:

$\frac{y_{i+1}-y_i}{\tau}=f(y_{i+1})$

## Implementation

The implicit euler method requires solving a nonlinear system at each time step
since the function evaluation is at the unknown future point $y_{i+1}$.
This is done using a Newton solver.

```
m_yold->set(y);
m_tau->set(tau);
NewtonSolver(m_equ, y);
```

## Results for Mass-Spring System

Below is a phase plot for a mass-spring system solved using the Implicit Euler method.

![Mass Spring Phase Plot Implicit Euler Method](../../pictures/MassSpringPhasePlot_implicit.png)

![Mass Spring System Time Evolution Implicit Euler Method](../../pictures/MassSpringSystemTimeEvolution_implicit.png)

