import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import csv

# Constants
R = 8314.0
P = 101325.0 # Pa
PG = 1 # atm
T = 1033  # System temperature (K)
RG = 82.06  # cm³ atm/ mol K

# Functions
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
    Q = (RG * T / PG) * N_total

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

# Initial conditions
B0 = np.array([1.0, 0.0, 0.0, 0.0])

# Volume range
V_start = 0.0
V_end = 1500.0

# Solver configuration
sol = solve_ivp(
    fun=decomposicao_benezo,
    t_span=(V_start, V_end),
    y0=B0,
    method='Radau',
    t_eval=np.linspace(V_start, V_end, 1000),
    rtol=1e-8,
    atol=1e-8
)

#    # Save results to CSV
# with open("output_stiff_benzeno.csv", "w", newline="") as csvfile:
#     writer = csv.writer(csvfile)
#     #writer.writerow(["Volume"] + [f"N{i+1}" for i in range(len(N0))])
#     for v, n in zip(sol.t, sol.y.T):
#         writer.writerow([v] + n.tolist())

# # Display a preview of results
# print("Volume (V), ", ", ".join([f"N{i+1}" for i in range(len(B0))]))
# for v, n in zip(sol.t[-10:], sol.y.T[-10:]):  # Display the first 10 rows
#     print(f"{v:.2f}, " + ", ".join([f"{ni:.5e}" for ni in n]))


# Plot results
plt.figure(figsize=(10, 6))
labels = ["C6H6", "H", "C12H10", "C18H14"]
for i in range(len(B0)):
    plt.plot(sol.t, sol.y[i], label=labels[i])

plt.xlabel("Volume (cm³)")
plt.ylabel("Molar Flow (mol/s)")
plt.title("Decomposição do Benzeno - Perfis de Concentração")
plt.legend()
plt.grid()
plt.show()
