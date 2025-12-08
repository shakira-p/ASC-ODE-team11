# Runge Kutta Methods

## Introduction

The Runge Kutta methods are a family of iterative methods for solving ordinary differential equations (ODEs). 
They provide a powerful framework for achieving higher-order accuracy in time-stepping schemes. 
This document describes both explicit and implicit Runge Kutta methods implemented in the ASC-ODE package.

## Mathematical Overview

Given an ODE $\frac{dy}{dt} = f(t, y)$, the Runge Kutta methods compute the solution at time $t_{n+1} = t_n + \tau$ using a series of intermediate stages.

Runge–Kutta methods solves the initial value problem
$y'(t) = f(t, y(t)), \quad y(t_0) = y_0$
by approximating
$y(t_{n+1}) = y(t_n) + \int_{t_n}^{t_{n+1}} f(s, y(s)) \, ds$
with a weighted sum of a few function evaluations inside the step. [web:2]

In an $s$‑stage RK method, one computes intermediate slopes
$k_j \approx f(t_n + c_j h,\, y_n + h \sum_{\ell} a_{j\ell} k_\ell)$
and updates
$y_{n+1} = y_n + h \sum_{j=1}^s b_j k_j$,  
where $h$ is the step size and $a_{j\ell}, b_j, c_j$ are method‑defining coefficients. [web:2]

These coefficients are usually arranged in a Butcher tableau,
and different choices give different methods, such as the explicit midpoint method, 
the classical RK4 method, or implicit Gauss/Radau schemes with better stability.

## Implementation

Both explicit and implicit Runge Kutta methods are implemented in the ASC-ODE package.

### Explicit Runge Kutta Method

The explicit Runge Kutta method is implemented using the Butcher tableau representation. 
The method computes intermediate stages based on the coefficients defined in the tableau and 
combines them to update the solution. It can be see as an extension of the explicit Euler method to higher-order accuracy.

#### Constructor Parameters

- `rhs`: A `NonlinearFunction` representing $f(t, y)$
- `A`: Coefficient matrix from the Butcher tableau
- `b`: Weight vector from the Butcher tableau
- `c`: Node vector from the Butcher tableau

#### DoStep Parameters

- `tau`: Time step size $\tau$
- `y`: Current solution vector

#### Key Features

- Computes the solution using a fixed sequence of explicit stages where each stage depends only on previously computed stages, so no nonlinear or linear systems must be solved.
- Simple to implement and flexible: different orders and error properties are obtained by choosing Butcher‑tableau coefficients, making it easy to build higher‑order schemes.
- Conditionally stable and best suited for non‑stiff problems — time step size must satisfy stability constraints (e.g., CFL‑type limits) for reliable results
### Implicit Runge Kutta Method

The implicit Runge Kutta method is also implemented using the Butcher tableau representation.
It gives optimal accuracy because it is based on the Gaussian quadrature points.

#### Constructor Parameters

- `rhs`: A `NonlinearFunction` representing $f(t, y)$
- `A`: Coefficient matrix from the Butcher tableau
- `b`: Weight vector from the Butcher tableau
- `c`: Node vector from the Butcher tableau


#### DoStep Parameters

- `tau`: Time step size $\tau$
- `y`: Current solution vector

#### Key Features

- Builds the nonlinear system for the concatenated stage vector $k = [k_0;...;k_{s-1}]: F(k) = k - f( y_n + tau * (A ⊗ I) * k ) = 0$
- Uses NewtonSolver(m_equ, m_k) to solve for k, then updates: $y_{n+1} = y_n + \tau * sum_j b_j * k_j$
- Stores m_k and m_y as contiguous vectors: stage j occupies indices $[j*n, (j+1)*n)$.

### Helper Functions

The package includes helper functions to generate common Butcher tableaus for both explicit and 
implicit Runge Kutta methods.

`Gauss2a,Gauss2b, Gauss2c, Gauss3c`
- Purpose: Predefined Butcher coefficients / nodes for common Gauss–Legendre RK schemes (2-stage, 3-stage nodes).

`GaussLegendre(VectorView<> x, VectorView<> w)`
- Purpose: Compute n-point Gauss–Legendre quadrature nodes x and weights w on [0,1].
- Inputs/Outputs: x and w are mutable views of length n; filled in-place.

`GaussJacobi(VectorView<> x, VectorView<> w, double alf, double bet)`
- Purpose: Compute n-point Gauss–Jacobi nodes/weights for parameters α,β.
- Inputs/Outputs: x, w views; alf, bet parameters.

`ComputeABfromC(const Vector<> &c)`
- Purpose: Given collocation nodes c, compute Butcher matrix A and weight vector b so the collocation RK satisfies moment conditions.
- Inputs/Outputs: c (vector) → returns (A,b) (Matrix, Vector).
- Algorithm: builds Vandermonde-like matrix M of powers of c, inverts M, uses integrals 1/(i+1) to get b, and similarly for rows of A.

`GaussRadau(VectorView<> x, VectorView<> w)`
- Purpose: Compute Gauss–Radau nodes/weights on [0,1] with the last node fixed at 1.
- Implementation: Calls GaussJacobi for n-1 inner nodes (α=1,β=0), rescales to [0,1], sets x[n-1]=1, and computes last weight to sum to 1.

## Results for Mass-Spring System

Below is a phase plot for a mass-spring system solved using both Explicit and Implicit Runge Kutta methods with the different helper functions to generate Butcher tableaus.

First the system phase plot and time evolution using the Explicit Runge Kutta method with step size h=0.1:

![Mass Spring Phase Plot Explicit Runge Kutta Method](../../pictures/MassSpringPhasePlot_rungekutta_explicit_100.png)
![Mass Spring System Time Evolution Explicit Runge Kutta Method](../../pictures/MassSpringSystemTimeEvolution_rungekutta_explicit_100.png)

Next the system phase plot and time evolution using the Implicit Runge Kutta method with step size h=0.1:

![Mass Spring Phase Plot Implicit Runge Kutta Method](../../pictures/MassSpringPhasePlot_rungekutta_implicit_100.png)
![Mass Spring System Time Evolution Implicit Runge Kutta Method](../../pictures/MassSpringSystemTimeEvolution_rungekutta_implicit_100.png)

TODO add rest of mass spring plots system radau, legrendre, gauss