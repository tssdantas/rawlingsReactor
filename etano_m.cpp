#include "mex.hpp"
#include "mexAdapter.hpp"
#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>

using namespace matlab::data;
using matlab::mex::ArgumentList;

// Constants
const double A1 = 1e14;
const double E1 = 217600.0;
const double A2 = 3e14;
const double E2 = 165300.0;
const double A3 = 3.4e12;
const double E3 = 28500.0;
const double A4 = 1e12;
const double E4 = 0.0;
const double A5 = 1e13;
const double E5 = 200800.0;
const double A6 = 1e12;
const double E6 = 0.0;
const double R = 8.314;
const double P = 101325.0;

// Function to calculate reaction rate constant
double calculate_k(double A, double E, double T) {
    return A * exp(-E / (R * T));
}

// Function to calculate concentration
double calculate_concentration(double N_species, double N_total, double T) {
    return (P / (R * T)) * (N_species / N_total);
}

// ODE system function
void decomposicao_etano(const std::vector<double>& N, std::vector<double>& dNdV, double V, double T) {
    double k1 = calculate_k(A1, E1, T);
    double k2 = calculate_k(A2, E2, T);
    double k3 = calculate_k(A3, E3, T);
    double k4 = calculate_k(A4, E4, T);
    double k5 = calculate_k(A5, E5, T);
    double k6 = calculate_k(A6, E6, T);

    double N_total = std::accumulate(N.begin(), N.end(), 0.0);

    double C_c2h6 = calculate_concentration(N[0], N_total, T);
    double C_no = calculate_concentration(N[1], N_total, T);
    double C_c2h5 = calculate_concentration(N[2], N_total, T);
    double C_hno = calculate_concentration(N[3], N_total, T);
    double C_h = calculate_concentration(N[4], N_total, T);

    double r1 = k1 * C_c2h6 * C_no;
    double r2 = k2 * C_c2h5;
    double r3 = k3 * C_h * C_c2h6;
    double r4 = k4 * C_h * C_no;
    double r5 = k5 * C_hno;
    double r6 = k6 * C_c2h5 * C_hno;

    dNdV[0] = -r1 - r3 + r6;
    dNdV[1] = -r1 - r4 + r5 + r6;
    dNdV[2] = r1 - r2 + r3 - r6;
    dNdV[3] = r1 + r4 - r5 - r6;
    dNdV[4] = r2 - r3 - r4 + r5;
    dNdV[5] = r2;
    dNdV[6] = r3;
}

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        checkArguments(inputs);

        const double T = 1050.0; // Fixed temperature
        const double V_start = 0.0;
        const double V_end = 0.0015;
        const size_t num_steps = 1000;
        const double dV = (V_end - V_start) / num_steps;

        std::vector<double> N = {6.62e-3, 3.48e-4, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::vector<std::vector<double>> result;

        for (size_t step = 0; step <= num_steps; ++step) {
            double V = V_start + step * dV;
            std::vector<double> dNdV(7, 0.0);

            decomposicao_etano(N, dNdV, V, T);

            // Euler integration
            for (size_t i = 0; i < N.size(); ++i) {
                N[i] += dV * dNdV[i];
            }

            result.push_back(N);
        }

        // Convert result to MATLAB array
        ArrayFactory factory;
        TypedArray<double> output = factory.createArray<double>({result.size(), 7});
        for (size_t i = 0; i < result.size(); ++i) {
            for (size_t j = 0; j < 7; ++j) {
                output[i][j] = result[i][j];
            }
        }

        outputs[0] = output;
    }

private:
    void checkArguments(const ArgumentList& inputs) {
        if (inputs.size() != 0) {
            throw std::invalid_argument("This function does not take any input arguments.");
        }
    }
};
