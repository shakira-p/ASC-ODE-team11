# Runge Kutta Methods

## Introduction

## Mathematical Overview
 

## Implementation

Both explicit and implicit Runge Kutta methods are implemented in the ASC-ODE package.

### Explicit Runge Kutta Method

#### Constructor Parameters

#### DoStep Parameters

#### Key Features

```

```

### Implicit Runge Kutta Method

#### Constructor Parameters

#### DoStep Parameters

#### Key Features

- Builds the nonlinear system for the concatenated stage vector $k = [k_0;...;k_{s-1}]: F(k) = k - f( y_n + tau * (A ⊗ I) * k ) = 0$
- Uses NewtonSolver(m_equ, m_k) to solve for k, then updates: $y_{n+1} = y_n + \tau * sum_j b_j * k_j$
- Stores m_k and m_y as contiguous vectors: stage j occupies indices $[j*n, (j+1)*n)$.

### Helper Functions

 

