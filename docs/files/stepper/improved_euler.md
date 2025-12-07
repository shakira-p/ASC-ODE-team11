# Improved Euler 

## Definition

$y^{\sim} = y_n + \frac{\tau}{2}f(y_n)$
$y_{n+1} = y_n + \tau f(y^{\sim})$

## Code 
In the improved euler method, we first evaluate the function at the current point y to get an estimate of the slope. We then use this slope to compute an intermediate point y_hat halfway through the time step. Finally,
we evaluate the function at this intermediate point and use that value to update y for the full time step. 
We continue this process for each time step until we reach the desired end time.

Here's how the evalutation in the timestepper is implemented in C++:

```
//evaluating f at current y
this->m_rhs->evaluate(y, m_vecf);

// intermediate point
m_y_hat = y;
m_y_hat += tau/2 * m_vecf;

// evaluating at y_hat
this->m_rhs->evaluate(m_y_hat, m_vecf);

// actual update: y += tau * f(y_hat)
y += tau * m_vecf;
```

## Results

Below is a phase plot for a mass-spring system solved using the Improved Euler method.

First the system phase plot and time evolution using the Improved Euler method with step size h=0.1:

![Mass Spring Phase Plot Improved Euler Method](/docs/plots/MassSpringPhasePlot_improved.png)

![Mass Spring System Time Evolution Improved Euler Method](/docs/plots/MassSpringSystemTimeEvolution_improved.png)
