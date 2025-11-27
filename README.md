# ASC-ODE
A package for solving ordinary differential equations

Read the [documentation](https://tuwien-asc.github.io/ASC-ODE/intro.html)

Find theory behind here: https://jschoeberl.github.io/IntroSC/ODEs/ODEs.html

# get the program running wuhu

do this in src folder:

```
mkdir build
cd build
cmake ..
make
```

this is what I needed to install:
```

sudo apt install python3.12-ven
sudo  apt install python3-pybind11

```
and I also needed to install cmake, but I think only one of them worked. Probaly the second one:
```
sudo apt install cmake  # version 4.1.2
sudo apt  install cmake  # version 3.27.8-1build
```

then itshould run with the first block of code

That worked thanks! :D


# get documentation running on Windows WSL wuhu

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