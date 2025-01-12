import numpy as np
from scipy.integrate import solve_ivp
import csv

# Constants
mu = 1000.0
RG1 = 8.314  # J/mol K
RG2 = 82.06  # cm³ atm/ mol K
T = 1050.0  # K
P = 1.0  # atm
Qf = 600.0  # cm³/s

A11 = 1e14
A12 = 1e12
A2 = 3e14
A3 = 3.4e12
A41 = 1e12
A42 = 1e13

E11 = 217600.0
E12 = 0.0
E2 = 165300.0
E3 = 28500.0
E41 = 0.0
E42 = 200800.0

# Functions
def calculate_k(A, E, T):
    return A * np.exp(-E / (RG1 * T))

def calculate_concentration(N_species, Q):
    return N_species / Q

# ODE system
def decomposicao_etano(V, N):
    # Calculate reaction rates
    k11 = calculate_k(A11, E11, T)
    k12 = calculate_k(A12, E12, T)
    k2 = calculate_k(A2, E2, T)
    k3 = calculate_k(A3, E3, T)
    k41 = calculate_k(A41, E41, T)
    k42 = calculate_k(A42, E42, T)
    
    # Total molar flow
    N_total = np.sum(N)
    Q = (RG2 * T / P) * N_total
    
    # Concentrations
    C_c2h6 = calculate_concentration(N[0], Q)
    C_c2h5_p = calculate_concentration(N[1], Q)
    C_c2h4 = calculate_concentration(N[2], Q)
    C_h = calculate_concentration(N[3], Q)
    C_h2 = calculate_concentration(N[4], Q)
    C_NO = calculate_concentration(N[5], Q)
    C_HNO = calculate_concentration(N[6], Q)
    
    # Reaction rates
    r1 = (k11 * C_c2h6 * C_NO) - (k12 * C_c2h5_p * C_HNO)
    r2 = k2 * C_c2h5_p
    r3 = k3 * C_h * C_c2h6
    r4 = (k41 * C_h * C_NO) - (k42 * C_HNO)
    
    # Differential equations
    dNdV = np.zeros_like(N)
    dNdV[0] = -r1 - r3
    dNdV[1] = r1 - r2 + r3
    dNdV[2] = r2
    dNdV[3] = r2 - r3 - r4
    dNdV[4] = r3
    dNdV[5] = -r1 - r4
    dNdV[6] = r1 + r4
    
    return dNdV

# Initial conditions
N0 = np.zeros(7)
N0[0] = (0.95 * Qf * P) / (RG2 * T)  # Initial concentration of C2H6
N0[5] = (0.05 * Qf * P) / (RG2 * T)  # Initial concentration of NO

# Volume range
V_start = 0.0
V_end = 1500.0

# Solver configuration
sol = solve_ivp(
    fun=decomposicao_etano,     # ODE function
    t_span=(V_start, V_end),   # Integration interval
    y0=N0,                     # Initial conditions
    method='BDF',              # Stiff solver
    t_eval=np.linspace(V_start, V_end, 1000),  # Evaluation points
    rtol=1e-6,                 # Relative tolerance
    atol=1e-8                  # Absolute tolerance
)

# Save results to CSV
with open("output_stiff.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Volume"] + [f"N{i+1}" for i in range(len(N0))])
    for v, n in zip(sol.t, sol.y.T):
        writer.writerow([v] + n.tolist())

# Display a preview of results
print("Volume (V), ", ", ".join([f"N{i+1}" for i in range(len(N0))]))
for v, n in zip(sol.t[:10], sol.y.T[:10]):  # Display the first 10 rows
    print(f"{v:.2f}, " + ", ".join([f"{ni:.5e}" for ni in n]))
