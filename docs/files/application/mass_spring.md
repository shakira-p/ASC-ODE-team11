# Mass Spring System

This example demonstrates the application of the ASC-ODE package to simulate
a mass-spring system. The mass-spring system is a classic example in physics
and engineering, characterized by oscillatory motion due to the restoring force
of the spring.

## Explicit Euler

Below is a phase plot for a mass-spring system solved using the Explicit Euler method.

![Mass Spring Phase Plot Explicit Euler Method](../../pictures/MassSpringPhasePlot_explict.png)

![Mass Spring System Time Evolution Explicit Euler Method](../../pictures/MassSpringSystemTimeEvolution_explict.png)

Here we can see that the Explicit Euler method does not conserve energy,
leading to a spiraling outwards trajectory in the phase plot.

## Implicit Euler

Below is a phase plot for a mass-spring system solved using the Implicit Euler method.

![Mass Spring Phase Plot Implicit Euler Method](../../pictures/MassSpringPhasePlot_implicit.png)

![Mass Spring System Time Evolution Implicit Euler Method](../../pictures/MassSpringSystemTimeEvolution_implicit.png)

Here we observe that the Implicit Euler method also fails to conserve energy, 
resulting in a spiraling inwards trajectory in the phase plot.

## Improved Euler

Below is a phase plot for a mass-spring system solved using the Improved Euler method.

First the system phase plot and time evolution using the Improved Euler method with step size h=0.1:

![Mass Spring Phase Plot Improved Euler Method](../../pictures/MassSpringPhasePlot_improved.png)

![Mass Spring System Time Evolution Improved Euler Method](../../pictures/MassSpringSystemTimeEvolution_improved.png)

Here we can see that the Improved Euler method provides a better approximation of the system's dynamics,
with a more circular trajectory in the phase plot, indicating improved energy conservation compared to the Explicit and Implicit Euler methods.

## Crank-Nicolson

Below is a phase plot for a mass-spring system solved using the Crank-Nicolson method.  

First the system phase plot and time evolution using the Improved Euler method with step size h=0.1:

![Mass Spring Phase Plot Improved Euler Method](../../pictures/MassSpringPhasePlot_crank.png)

![Mass Spring System Time Evolution Improved Euler Method](../../pictures/MassSpringSystemTimeEvolution_crank.png)

Here we can see that the Crank-Nicolson method provides a very good approximation of the system's dynamics,
with a nearly circular trajectory in the phase plot, indicating excellent energy conservation.


## Runge-Kutta

### Explicit Runge Kutta Method

Below is a phase plot for a mass-spring system solved using both Explicit and Implicit Runge Kutta methods with the different helper functions to generate Butcher tableaus.

First the system phase plot and time evolution using the Explicit Runge Kutta method with step size h=0.1:

![Mass Spring Phase Plot Explicit Runge Kutta Method](../../pictures/MassSpringPhasePlot_rungekutta_explicit_100.png)
![Mass Spring System Time Evolution Explicit Runge Kutta Method](../../pictures/MassSpringSystemTimeEvolution_rungekutta_explicit_100.png)


### Implicit Runge Kutta Method

Next the system phase plot and time evolution using the Implicit Runge Kutta method with step size h=0.1:

![Mass Spring Phase Plot Implicit Runge Kutta Method](../../pictures/MassSpringPhasePlot_rungekutta_implicit_100.png)
![Mass Spring System Time Evolution Implicit Runge Kutta Method](../../pictures/MassSpringSystemTimeEvolution_rungekutta_implicit_100.png)

#### Implicit Runge Kutta Method with Gauss-Legendre Butcher Tableau

Next the system phase plot and time evolution using the Implicit Runge Kutta method with Gauss-Legendre Butcher Tableau and step size h=0.1:

![Mass Spring Phase Plot Implicit Runge Kutta Method Gauss-Legendre](../../pictures/MassSpringPhasePlot_rungekutta_gausslegrendre_100.png)
![Mass Spring System Time Evolution Implicit Runge Kutta Method Gauss-Legendre](../../pictures/MassSpringSystemTimeEvolution_rungekutta_gausslegrendre_100.png)

## Comparison of different methods and step sizes

### Step Size T = 100

![Mass Spring Phase Plot Step Size 100](../../pictures//MassSpringPhasePlot_100.png)
![Mass Spring System Time Evolution Step Size 100](../../pictures//MassSpringSystemTimeEvolution_100.png)

### Step Size T = 150

![Mass Spring Phase Plot Step Size 100](../../pictures//MassSpringPhasePlot_150.png)
![Mass Spring System Time Evolution Step Size 100](../../pictures//MassSpringSystemTimeEvolution_150.png)
### Step Size T = 200

![Mass Spring Phase Plot Step Size 100](../../pictures//MassSpringPhasePlot_200.png)
![Mass Spring System Time Evolution Step Size 100](../../pictures//MassSpringSystemTimeEvolution_200.png)