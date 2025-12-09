# Introduction to ASC-ODE

ASC-ODE is a C++ library for solving ordinary differential equations (ODEs).
The equation is defined by the right hand side function.
ASC-ODE provides various time-steppers which may be used for odes with right hand sides
given by a function object.

Find theory behind here: https://jschoeberl.github.io/IntroSC/ODEs/ODEs.html

## How to run the libarys

Make sure you have the necessary dependencies installed, such as Python 3.12 and pybind11.
You can install them using:

```sh
sudo apt install python3.12-ven
sudo apt install python3-pybind11
``` 

To build the library, navigate to the `src` folder and execute the following commands:

```sh
mkdir build
cd build
cmake .. 
make
```

After the build is complete, you can run the demonstration scripts. There should be an executable named `test_odes` in the `build` directory.

To run the demonstration scripts, execute the following command from the `build` directory:
```sh
./test_odes {step_size_of_your_choice}
```

This will give you a text output of the simulation results. This needs to be copied to the `demo/data` folder for further processing.

To run the mass-spring system simulation and generate plots with python, navigate in the demo folder and use the following command:

```sh
python3 ./plot_mass_spring.py
```

### Change simulation method

To change the simulation method, you can modify the `test_ode.cpp` file located in the `src` folder.
There are different time-stepping methods available, such as Explicit Euler, Implicit Euler, Improved Euler, Crank-Nicolson, and Runge-Kutta methods.
You can uncomment the desired method and comment out the others to switch between them.
For example, to use the Explicit Euler method, you would modify the code as follows:

```cpp
// Uncomment the Explicit Euler method
ExplicitEuler stepper(rhs)
// Comment out the other methods
// ImplicitEuler stepper(rhs);
// ImprovedEuler stepper(rhs);
// CrankNicolson stepper(rhs);

std::ofstream outfile ("../rungekutta_explicit.txt");

```

After making the changes, rebuild the library and run the demonstration scripts again to see the results with the new method.
There can only be one method active at a time. So make sure to comment out the other methods when switching.
If you want to try multiple methods, you can run the program multiple times with different methods uncommented each time.
Make sure to change the output file name accordingly and name it with your method and step size to avoid overwriting previous results.