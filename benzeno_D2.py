import numpy as np
from scipy.integrate import solve_ivp
import csv

# Constants
mu = 1000.0
A1 = 1e14
E1 = 217600.0
A2 = 3e14
E2 = 165300.0
A3 = 3.4e12
E3 = 28500.0
A4 = 1e12
E4 = 0.0
A5 = 1e13
E5 = 200800.0
A6 = 1e12
E6 = 0.0
R = 8314.0
P = 101325.0
T = 1033  # System temperature (K)

# Functions
def calculate_k(A, E, T):
    return A * np.exp(-E / (R * T))

def calculate_concentration(N_species, N_total, T):
    return (P / (R * T)) * (N_species / N_total)

def decomposicao_benezo(V, N):
    """ODE system for benzene decomposition."""
    _k1 = 7.0e5  # L/mol.h
    K1 = 0.31
    _k2 = 4e5  # L/mol.h
    K2 = 0.48

    # Total molar flow
    N_total = np.sum(N)

    # Concentrations
    C_c6h6 = calculate_concentration(N[0], N_total, T)
    C_h = calculate_concentration(N[1], N_total, T)
    C_c12h10 = calculate_concentration(N[2], N_total, T)
    C_c18h14 = calculate_concentration(N[3], N_total, T)

    # Reaction rates
    r1 = _k1 * (C_c6h6**2 - (C_c12h10 * C_h / K1))
    r2 = _k2 * (C_c6h6 * C_c12h10 - (C_c18h14 * C_h / K2))

    # Differential equations
    dNdV = np.zeros_like(N)
    dNdV[0] = -(2 * r1) - r2    # dNc2h6/dV
    dNdV[1] = r1 + r2           # dNno/dV
    dNdV[2] = r1 - r2           # dNc2h5/dV
    dNdV[3] = r2                # dNhno/dV

    return dNdV

# Observer
def my_observer(V, N, file):
    """Observer function to log results."""
    file.writerow([V] + N.tolist())
    #print(f"V: {V:.2f}\t" + "\t".join([f"{x:.5e}" for x in N]))

# Main integration routine
if __name__ == "__main__":
    # Initial conditions
    
    B0 = np.array([1.0, 0.0, 0.0, 0.0])

    # Volume range
    V_start = 0.0
    V_end = 1500.0
    V_step = (V_end - V_start) / 1e4

    # Solver configuration
    sol = solve_ivp(
        fun=decomposicao_benezo,     # ODE function
        t_span=(V_start, V_end),   # Integration interval
        y0=B0,                     # Initial conditions
        method='LSODA',              # Stiff solver BDF, Radau, LSODA 
        t_eval=np.linspace(V_start, V_end, 1000),  # Evaluation points
        rtol=1e-6,                 # Relative tolerance
        atol=1e-8                  # Absolute tolerance
    )

    # Save results to CSV
    with open("output_stiff_benzeno.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        #writer.writerow(["Volume"] + [f"N{i+1}" for i in range(len(N0))])
        for v, n in zip(sol.t, sol.y.T):
            writer.writerow([v] + n.tolist())

    # Display a preview of results
    print("Volume (V), ", ", ".join([f"N{i+1}" for i in range(len(B0))]))
    for v, n in zip(sol.t[-10:], sol.y.T[-10:]):  # Display the first 10 rows
        print(f"{v:.2f}, " + ", ".join([f"{ni:.5e}" for ni in n]))
