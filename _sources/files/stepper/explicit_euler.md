# Explicit Euler

## Definition

$y_{i+1} = y_i + \tau f(y_i)$ $\space  0 \leq i\leq n$ with  $\tau = \frac{T}{n}$

## Code

The explicit euler method is the simplest of the ODE solvers implemented in the ASC-ODE package. It updates the solution by taking a step in the direction of the derivative at the current point.

```
this->m_rhs->evaluate(y, m_vecf);
y += tau * m_vecf;
```

## Results

Below is a phase plot for a mass-spring system solved using the Explicit Euler method.

![Mass Spring Phase Plot Explicit Euler Method](/docs/plots/MassSpringPhasePlot_explict.png)

![Mass Spring System Time Evolution Explicit Euler Method](/docs/plots/MassSpringSystemTimeEvolution_explict.png)

Here we can see that the Explicit Euler method does not conserve energy, leading to a spiraling outwards trajectory in the phase plot.