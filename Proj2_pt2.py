import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt

# Lagrange basis functions
def lagrange_basis(xi, nodes, i):
    basis = 1.0
    for j in range(len(nodes)):
        if j != i:
            basis *= (xi - nodes[j]) / (nodes[i] - nodes[j])
    return basis

def lagrange_pbasis(xi, nodes, i):
    basis_prime = 0.0
    for j in range(len(nodes)):
        if j != i:
            basis_prime += 1 / (nodes[i] - nodes[j]) * lagrange_basis(xi, nodes[:j] + nodes[j+1:], i)
    return basis_prime

# Assemble global matrices
def assemble_global_matrices(N, nodes, f_func, xl, xr, dt):
    Ne = N - 1
    K_Gmtrx = np.zeros((N, N))  # Global stiffness matrix
    M_Gmtrx = np.zeros((N, N))  # Global mass matrix
    F_Gvec = np.zeros(N)        # Global load vector

    for k in range(Ne):
        nodes_local = nodes[k:k+2]  # Nodes of the local element
        det_J = (xr - xl) / 2       # Jacobian determinant for mapping from local to global

        K_Lmtrx = np.zeros((2, 2))  # Local stiffness matrix
        M_Lmtrx = np.zeros((2, 2))  # Local mass matrix
        F_Lvec = np.zeros(2)        # Local load vector

        for l in range(2):
            for m in range(2):
                # Elemental stiffness matrix
                integrand_K = lambda xi: lagrange_pbasis(xi, nodes_local, l) * lagrange_pbasis(xi, nodes_local, m) * det_J
                K_Lmtrx[l, m] = np.trapz([integrand_K(xi) for xi in nodes_local], nodes_local)  # Numerical integration using trapezoidal rule

                # Elemental mass matrix
                integrand_M = lambda xi: lagrange_basis(xi, nodes_local, l) * lagrange_basis(xi, nodes_local, m) * det_J
                M_Lmtrx[l, m] = np.trapz([integrand_M(xi) for xi in nodes_local], nodes_local)  # Numerical integration using trapezoidal rule

            # Elemental load vector
            integrand_F = lambda xi: f_func(xi, dt) * lagrange_basis(xi, nodes_local, l) * det_J
            F_Lvec[l] = np.trapz([integrand_F(xi) for xi in nodes_local], nodes_local)          # Numerical integration using trapezoidal rule

        # Assemble local to global matrices
        for l in range(2):
            global_node1 = k + l
            F_Gvec[global_node1] += F_Lvec[l]

            for m in range(2):
                global_node2 = k + m
                K_Gmtrx[global_node1, global_node2] += K_Lmtrx[l, m]
                M_Gmtrx[global_node1, global_node2] += M_Lmtrx[l, m]

    return F_Gvec, K_Gmtrx, M_Gmtrx

# Solve heat equation using backward Euler
def HeatEq_Bkwd(K, M, dt, u):
    A = M + dt * K
    inv_A = np.linalg.inv(A)
    u = np.dot(inv_A, np.dot(M, u))
    return u

# User Input:
def user_func(xi, dt):
    return (np.pi**2 - 1) * np.exp(-dt) * np.sin(np.pi * xi)

# Number of nodes
N = 11
# Domain boundaries
xl = 0
xr = 1
# Initial and final time
T0 = 0
Tf = 0.005
# Time step
dt = 0.001
# Initialize nodes
nodes = np.linspace(xl, xr, N)
# Set initial condition
u0 = np.sin(np.pi * nodes)

# Spatial step size
dx = (xr - xl) / (N - 1)

# Time-stepping with backward Euler
nt = int((Tf - T0) / dt)
u_curr = u0
for n in range(1, nt + 1):
    ctime = T0 + n * dt

    # Assemble global matrices
    F_Gvec, K_Gmtrx, M_Gmtrx = assemble_global_matrices(N, nodes, user_func, xl, xr, ctime)

    # Solve the heat equation using backward Euler
    u_curr = HeatEq_Bkwd(K_Gmtrx, M_Gmtrx, dt, u_curr)

# Plot the results
plt.plot(nodes, u_curr, label='Backward Euler')
plt.plot(nodes, np.exp(-Tf) * np.sin(np.pi * nodes), label='Analytical Solution', linestyle='--')
plt.title('Solution with Backward Euler at Final Time')
plt.xlabel('X')
plt.ylabel('Displacement')
plt.legend()
plt.show()
