# COE352_Proj2
This repository is for my COE352 project, which involves solving the heat equation using the finite element method with forward and backward Euler time-stepping.  

# Heat Transfer Problem

The heat transfer problem is described by the equation: `uₜ - u_{ₓₓ} = f(x, t), (x, t) ∈ (0,1) × (0,1)`


The analytic solution to this problem is given by: ` u(x, t) = e^{-t}sin(πx) `

## Table Of Contents
- [Proj2_Handwritten](#proj2_handwritten)
  - [Description](#description)
- [Forward Euler Code](#proj1-pt1py)
  - [Description](#description-1)
  - [Usage](#usage)
  - [Functions](#functions)
  - [Example Usage](#example-usage)
  - [Free Response Answers & Plots](#free-response-answers--plots)
- [Backward Euler Code](#Proj1_pt2.py)
  - [Description](#description-2)
  - [Usage](#usage-1)
  - [Functions](#functions-1)
  - [Example Usage](#example-usage-1)
  - [Free Response Answers & Plots](#free-response-answers--plots-1)

## Proj2_Handwritten

### Description  
This image is the derived weak form of the above Heat Equation.

![Hand-written Weak Form](https://github.com/kjmadera24/PictureFiles/blob/aceba3751ed8050617254bec18c52b916b486a2e/File_000%20(4).png)

## Proj2_pt1.py

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

### Free Response Answers & Plots  

1. **Plotting the Results:** * from T final going from 0.24 to 0.00024. *
<p float="left">
  <img src="https://github.com/kjmadera24/PictureFiles/blob/aceba3751ed8050617254bec18c52b916b486a2e/Pt1%20Tf_0.24.png" width="400"/>
  <img src="https://github.com/kjmadera24/PictureFiles/blob/aceba3751ed8050617254bec18c52b916b486a2e/Pt1%20Tf_0.024.png" width="400"/>
  <img src="https://github.com/kjmadera24/PictureFiles/blob/aceba3751ed8050617254bec18c52b916b486a2e/Pt1%20Tf_0.0024.png" width="400"/>
  <img src="https://github.com/kjmadera24/PictureFiles/blob/aceba3751ed8050617254bec18c52b916b486a2e/Pt1%20Tf_0.00024.png" width="400"/>
</p>

2. **Question 1:** Increase the time-step until you find the instability, what dt does this occur at?
*The staibility increases as dt increases and decreases as dt decreases, with this code I tried from 1/240 up to 1/24000000 for dt and I still didnt see much difference. I am sure with a stronger laptop I would eventually see the instability more.*  
<p float="left">
  <img src="https://github.com/kjmadera24/PictureFiles/blob/aceba3751ed8050617254bec18c52b916b486a2e/Pt1%20dt_2400.png" width="400" width="400"/>
  <img src="https://github.com/kjmadera24/PictureFiles/blob/aceba3751ed8050617254bec18c52b916b486a2e/Pt1%20dt_2400.png" width="400"/>
  <img src="https://github.com/kjmadera24/PictureFiles/blob/aceba3751ed8050617254bec18c52b916b486a2e/Pt1%20dt_24000.png" width="400"/>
  <img src="https://github.com/kjmadera24/PictureFiles/blob/aceba3751ed8050617254bec18c52b916b486a2e/Pt1%20dt_24000000.png" width="400"/>
</p>

3. **Question 2:** How does the solution change as N decreases?
*The plot gets more round the higher the N count. A way to think about it is that a square has 4 sides and is blocky, but a dodecahedron has 12 sides and looks more round.*  

<p float="left">
  <img src="https://github.com/kjmadera24/PictureFiles/blob/aceba3751ed8050617254bec18c52b916b486a2e/Pt1%20N_4.png" width="400"/>
  <img src="https://github.com/kjmadera24/PictureFiles/blob/aceba3751ed8050617254bec18c52b916b486a2e/Pt1%20N_8.png" width="400"/>
  <img src="https://github.com/kjmadera24/PictureFiles/blob/aceba3751ed8050617254bec18c52b916b486a2e/Pt1%20N_11.png" width="400"/>
  <img src="https://github.com/kjmadera24/PictureFiles/blob/aceba3751ed8050617254bec18c52b916b486a2e/Pt1%20N_24.png" width="400"/>
</p>

## Proj2_pt2.py

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

### Free Response Answers & Plots  

**Plotting the Results:** What happens as the time-step is equal to or greater than the spatial step size?  

*When you make the time step equal to or greater than the spacial time step the plotting gets more accurate. This is opposite to the Forward Euler, in these plots you can see how stable the plots get going from being a dt of 0.001 to 10.*

<p float="left">
  <img src="https://github.com/kjmadera24/PictureFiles/blob/01b1bd6b463bb27569a5d9468d1bc4f56a162aa8/Pt2%20dt_0.001.png" width="400"/>
  <img src="https://github.com/kjmadera24/PictureFiles/blob/01b1bd6b463bb27569a5d9468d1bc4f56a162aa8/Pt2%20dt_0.1.png" width="400"/>
  <img src="https://github.com/kjmadera24/PictureFiles/blob/01b1bd6b463bb27569a5d9468d1bc4f56a162aa8/Pt2%20dt_1.png" width="400"/>
  <img src="https://github.com/kjmadera24/PictureFiles/blob/01b1bd6b463bb27569a5d9468d1bc4f56a162aa8/Pt2%20dt_10.png" width="400"/>
</p>

