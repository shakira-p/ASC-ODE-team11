# Implicit Euler

## Definition

For every timestep for the value $y^{n+1}$ we need to solve:

$\frac{y_{i+1}-y_i}{\tau}=f(y_{i+1})$

## Code

The implicit euler method requires solving a nonlinear system at each time step
since the function evaluation is at the unknown future point $y_{i+1}$.
This is done using a Newton solver.

```
m_yold->set(y);
m_tau->set(tau);
NewtonSolver(m_equ, y);
```

## Results

Below is a phase plot for a mass-spring system solved using the Implicit Euler method.

![Mass Spring Phase Plot Implicit Euler Method](/docs/plots/MassSpringPhasePlot_implicit.png)

![Mass Spring System Time Evolution Implicit Euler Method](/docs/plots/MassSpringSystemTimeEvolution_implicit.png)

