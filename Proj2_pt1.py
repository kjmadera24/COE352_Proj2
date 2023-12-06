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

# Elemental matrices
def elmtl_mtrxs(N, nodes, f_func, xl, xr, dt):
    Ne = N - 1
    K_Gmtrx = np.zeros((N, N))  # Global stiffness matrix
    M_Gmtrx = np.zeros((N, N))  # Global mass matrix
    F_Gvec = np.zeros(N)        # Global load vector

    for k in range(Ne):
        nodes_local = nodes[k:k+2]  # Nodes of the local element
        det_J = (xr - xl) / 2       # Jacobian determinant for mapping from local to global coordinates

        K_Lmtrx = np.zeros((2, 2))  # Local stiffness matrix
        M_Lmtrx = np.zeros((2, 2))  # Local mass matrix
        F_Lvec = np.zeros(2)       # Local load vector

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
            F_Lvec[l] = np.trapz([integrand_F(xi) for xi in nodes_local], nodes_local)  # Numerical integration using trapezoidal rule

        # Assemble local matrices to global matrices
        for l in range(2):
            global_node1 = k + l
            F_Gvec[global_node1] += F_Lvec[l]

            for m in range(2):
                global_node2 = k + m
                K_Gmtrx[global_node1, global_node2] += K_Lmtrx[l, m]
                M_Gmtrx[global_node1, global_node2] += M_Lmtrx[l, m]

    return F_Gvec, K_Gmtrx, M_Gmtrx

# Solve heat equation using forward Euler
def HeatEq_Frwd(F, K, M, dt, u):
    inv_M = np.linalg.inv(M)
    u = u - dt * np.dot(inv_M, np.dot(K, u)) + dt * np.dot(inv_M, F)
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
Tf = 0.0024
# Time step
dt = 1/551
# Initialize nodes
nodes = np.linspace(xl, xr, N)
# Set initial condition
u0 = np.sin(np.pi * nodes)
# Time-stepping
nt = int((Tf - T0) / dt)
u_curr = u0

for n in range(1, nt + 1):
    ctime = T0 + n * dt
    # Elemental matrices
    F_Gvec, K_Gmtrx, M_Gmtrx = elmtl_mtrxs(N, nodes, user_func, xl, xr, ctime)
    # Solve the heat equation
    u_curr = HeatEq_Frwd(F_Gvec, K_Gmtrx, M_Gmtrx, dt, u_curr)

# Plot the results
plt.plot(nodes, u_curr, label='Final Displacement')
plt.plot(nodes, np.exp(-Tf) * np.sin(np.pi * nodes), label='Analytical Solution', linestyle='--')
plt.title('Solution with Forward Euler at Final Time')
plt.xlabel('X')
plt.ylabel('Displacement')
plt.legend()
plt.show()