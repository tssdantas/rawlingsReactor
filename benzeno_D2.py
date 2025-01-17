import numpy as np
from scipy.integrate import solve_ivp, ode
import csv
import matplotlib.pyplot as plt

# Constants
R = 8314.0
P = 101325.0 #Pa
PG = 1 #atm
T = 1033  # System temperature (K)
RG = 82.06  # cm³ atm/ mol K

# Functions
def calculate_concentration(N_species, N_total, T):
    return (P / (R * T)) * (N_species / N_total)

def calculate_concentration_Q(N_species, Q):
    return (N_species / Q)

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

    # C_c6h6 = calculate_concentration_Q(N[0], Q)
    # C_h = calculate_concentration_Q(N[1], Q)
    # C_c12h10 = calculate_concentration_Q(N[2], Q)
    # C_c18h14 = calculate_concentration_Q(N[3], Q)

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
def my_observer(V, N):
    """Observer function to log results."""
    #file.writerow([V] + N.tolist())
    observed_data.append((V, N.copy()))
    #print(f"V: {V:.2f}\t" + "\t".join([f"{x:.5e}" for x in N]))

# Main integration routine
if __name__ == "__main__":
    # Initial conditions
    
    B0 = np.array([1.0, 0.0, 0.0, 0.0])

    # Volume range
    V_start = 0.0
    V_end = 1500.0
    V_step = (V_end - V_start) / 1e10

    # Solver configuration
    # sol = solve_ivp(
    #     fun=decomposicao_benezo,     # ODE function
    #     t_span=(V_start, V_end),   # Integration interval
    #     y0=B0,                     # Initial conditions
    #     method='Radau',             # Stiff solver BDF, Radau, LSODA 
    #     t_eval=np.linspace(V_start, V_end, 10**6),  # Evaluation points
    #     rtol=1e-8,                 # Relative tolerance
    #     atol=1e-8                  # Absolute tolerance
    # )

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

    solver = ode(decomposicao_benezo)
    solver.set_integrator('dop853', atol=1e-9, rtol=1e-9)
    solver.set_initial_value(B0, 0.0)
    solver.set_solout(my_observer)

    observed_data = []

    # while solver.successful() and solver.t < V_end:
    #     solver.integrate(solver.t + V_step)

    solver.integrate(V_end)

    # Extract results for plotting
    V_vals = [data[0] for data in observed_data]
    N_vals = np.array([data[1] for data in observed_data])

    plt.figure(figsize=(10, 6))
    for i in range(N_vals.shape[1]):
        plt.plot(V_vals, N_vals[:, i], label=f"N[{i}]")
    plt.xlabel("Volume (V)")
    plt.ylabel("Molar Flow Rate (N)")
    plt.title("Molar Flow Rates vs. Volume")
    plt.legend()
    plt.grid()
    plt.show()

 

