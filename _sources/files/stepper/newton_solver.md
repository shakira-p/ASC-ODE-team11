# Newton Solver

## Introduction

The Newton solver is a numerical method for finding roots of nonlinear systems of equations. It iteratively refines an initial guess until the residual falls below a specified tolerance.

## Mathematical Overview

Given a nonlinear function $F: \mathbb{R}^n \to \mathbb{R}^m$, the Newton method finds $x^*$ such that $F(x^*) = 0$ using the iteration:

$$x_{k+1} = x_k - J_F(x_k)^{-1} F(x_k)$$

where $J_F(x_k)$ is the Jacobian matrix of $F$ evaluated at $x_k$.

The iteration continues until $\|F(x_k)\| < \text{tol}$ or the maximum number of steps is reached.

## Implementation

The `NewtonSolver` function implements this method with the following signature:

```cpp
void NewtonSolver(std::shared_ptr<NonlinearFunction> func, VectorView<double> x,
                  double tol = 1e-10, int maxsteps = 10,
                  std::function<void(int,double,VectorView<double>)> callback = nullptr)
```

## Parameters
- func: A NonlinearFunction object that provides evaluate() for <span>F(x)</span> and evaluateDeriv() for the Jacobian <span>J_F(x)</span>
- x: Initial guess vector (modified in-place to contain the solution)
- tol: Convergence tolerance for the residual norm (default: <span>10^{-10}</span>)
- maxsteps: Maximum number of Newton iterations (default: 10)
- callback: Optional function called at each iteration with (iteration, error, current_x)

## Behavior
The solver computes the Jacobian, inverts it, and updates the solution vector until convergence or reaching maxsteps. If convergence fails, a std::domain_error exception is thrown.
