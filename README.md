# ASC-ODE
A python package for solving ordinary differential equations, implemented in C++ and built with pybind11.

The following applications are already implemented:
* Mass Spring system
* Electric Networks

By using different the stepper functions:
* Implicit Euler
* Explicit Euler / Improved Euler
* Runge-Kutta methods
* Newton Solver
* Crank-Nicolson

And helper classes:
* AutoDiff for automatic exact differentiation
* MassSpringSystem with addable DistanceConstraints

Read the [documentation](https://shakira-p.github.io/ASC-ODE-team11/overview.html)!

Find theory behind here: https://jschoeberl.github.io/IntroSC/ODEs/ODEs.html

# Get the program running wuhu

do this in src folder:

```
mkdir build
cd build
cmake ..
make
```

this is what is needed to install:
```

sudo apt install python3.12-ven
sudo  apt install python3-pybind11

```
install cmake
```
sudo apt install cmake  # version 3.27.8-1build
```

# On Windows WSL, get the build running wuhu

```sh
# 1. Install pipx
sudo apt install pipx

# 2. Add it to your path (run this once)
pipx ensurepath

# 3. Install jupyter-book (creates the hidden venv automatically)
pipx install "jupyter-book<2.0"
pipx inject jupyter-book ghp-import

# 4. Open new terminal

# 5. Go to docs folder
cd docs/

# 6. Build
jupyter-book build .
```
