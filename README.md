# COE352_Proj2
This repository is for my COE352 project, which involves solving the heat equation using the finite element method with forward and backward Euler time-stepping.  

The heat transfer problem

\[ u_t - u_{xx} = f(x, t), \quad (x, t) \in (0, 1) \times (0, 1) \]

With initial and Dirichlet boundary conditions

\[ u(x, 0) = \sin(\pi x), \]

\[ u(0, t) = u(1, t) = 0 \]

and function

\[ f(x, t) = (\pi^2 - 1)e^{-t}\sin(\pi x) \]

The analytic solution to this problem is

\[ u(x, t) = e^{-t}\sin(\pi x) \]

- [Proj1_pt1.py](#proj1_pt1.py)
  - [Description](#description)
  - [Usage](#usage)
  - [Functions](#functions)
    - [`lagrange_basis(xi, nodes, i)`](#lagrange_basisxi-nodes-i)
    - [`lagrange_basis_prime(xi, nodes, i)`](#lagrange_basis_primexi-nodes-i)
    - [`elmtl_mtrxs(N, nodes, f_func, xl, xr, dt)`](#elmtl_mtrxsN-nodes-f_func-xl-xr-dt)
    - [`HeatEq_Frwd(F, K, M, dt, u)`](#heateq_frwdF-k-m-dt-u)
    - [`user_func(xi, dt)`](#user_funcxi-dt)
  - [Example Usage](#example-usage)
- [Proj1_pt2.py](#proj1-pt2py)
  - [Description](#description-1)
  - [Usage](#usage-1)
  - [Functions](#functions-1)
    - [`assemble_global_matrices(N, nodes, f_func, xl, xr, dt)`](#assemble_global_matricesN-nodes-f_func-xl-xr-dt)
    - [`HeatEq_Bkwd(K, M, dt, u)`](#heateq_bkwdK-m-dt-u)
    - [`f_example(xi, dt)`](#f_examplexi-dt)
  - [Example Usage](#example-usage-1)
  - [Requirements](#requirements)
  - [Testing](#testing)
  - [License](#license)
  - [Acknowledgments](#acknowledgments)


## Proj1_pt1.py

### Description
This script implements a forward Euler time-stepping solution for the 1D heat equation using finite element methods. It utilizes Lagrange basis functions and performs numerical integration to assemble global matrices for the stiffness, mass, and load. The solution is then obtained by solving the heat equation with forward Euler.

### Usage
1. Download the script.
2. Edit the `user_func` function to set the desired parameters such as the number of nodes, domain boundaries, time step, etc.
3. Run the script in a Python environment.

### Functions

#### `lagrange_basis(xi, nodes, i)`
Calculates the Lagrange basis function at a given point `xi`.

#### `lagrange_pbasis(xi, nodes, i)`
Calculates the derivative of the Lagrange basis function at a given point `xi`.

#### `elmtl_mtrxs(N, nodes, f_func, xl, xr, dt)`
Computes the elemental matrices for stiffness, mass, and load, and assembles them into global matrices.

#### `HeatEq_Frwd(F, K, M, dt, u)`
Solves the heat equation using forward Euler time-stepping.

#### `user_func(xi, dt)`
User-defined function for the source term in the heat equation.

### Example Usage

#### Set parameters
N = 11
xl, xr = 0, 1
T0, Tf = 0, 0
dt = 1/551

#### Initialize nodes and initial condition
nodes = np.linspace(xl, xr, N)
u0 = np.sin(np.pi * nodes)

#### Time-stepping
nt = int((Tf - T0) / dt)
u_curr = u0


## Proj1_pt2.py

### Description
This script implements a backward Euler time-stepping solution for the 1D heat equation using finite element methods. It utilizes Lagrange basis functions and performs numerical integration to assemble global matrices for the stiffness, mass, and load. The solution is then obtained by solving the heat equation with backward Euler.

### Usage
1. Download the script.
2. Edit the `user_func` function to set the desired parameters such as the number of nodes, domain boundaries, time step, etc.
3. Run the script in a Python environment.

### Functions

#### `lagrange_basis(xi, nodes, i)`
Calculates the Lagrange basis function at a given point `xi`.

#### `lagrange_basis_prime(xi, nodes, i)`
Calculates the derivative of the Lagrange basis function at a given point `xi`.

#### `assemble_global_matrices(N, nodes, f_func, xl, xr, dt)`
Computes the elemental matrices for stiffness, mass, and load, and assembles them into global matrices.

#### `HeatEq_Bkwd(K, M, dt, u)`
Solves the heat equation using backward Euler time-stepping.

### Example Usage

#### Set parameters
N = 11
xl, xr = 0, 1
T0, Tf = 0, 0
dt = 1/551

#### Initialize nodes and initial condition
nodes = np.linspace(xl, xr, N)
u0 = np.sin(np.pi * nodes)

#### Time-stepping with backward Euler
nt = int((Tf - T0) / dt)
u_curr = u0
