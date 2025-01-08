#include "mex.hpp"
#include "mexAdapter.hpp"
#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>

using namespace matlab::data;
using matlab::mex::ArgumentList;

// Constants
const double _k1 = 7.0e5; // L/mol.h
const double K1 = 0.31;
const double _k2 = 4.0e5; // L/mol.h
const double K2 = 0.48;
const double R = 8314.0;
const double P = 101325.0;

// Function to calculate concentration
double calculate_concentration(double N_species, double N_total, double T) {
    return (P / (R * T)) * (N_species / N_total);
}

// ODE system function
void decomposicao_benezo(const std::vector<double>& N, std::vector<double>& dNdV, double V, double T) {
    double N_total = std::accumulate(N.begin(), N.end(), 0.0);

    double C_c6h6 = calculate_concentration(N[0], N_total, T);
    double C_h = calculate_concentration(N[1], N_total, T);
    double C_c12h10 = calculate_concentration(N[2], N_total, T);
    double C_c18h14 = calculate_concentration(N[3], N_total, T);

    double r1 = _k1 * (std::pow(C_c6h6, 2) - ((C_c12h10 * C_h) / K1));
    double r2 = _k2 * ((C_c6h6 * C_c12h10) - ((C_c18h14 * C_h) / K2));

    dNdV[0] = -(2 * r1) - r2;    // dNc2h6/dV
    dNdV[1] = r1 + r2;           // dNno/dV
    dNdV[2] = r1 - r2;           // dNc2h5/dV
    dNdV[3] = r2;                // dNhno/dV
}

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        checkArguments(inputs);

        const double T = 1033.0; // Fixed temperature
        const double V_start = 0.0;
        const double V_end = 1500.0;
        const size_t num_steps = 1000000;
        const double dV = (V_end - V_start) / num_steps;

        std::vector<double> N = {1.0, 0.0, 0.0, 0.0}; // Initial conditions
        std::vector<std::vector<double>> result;

        for (size_t step = 0; step <= num_steps; ++step) {
            double V = V_start + step * dV;
            std::vector<double> dNdV(4, 0.0);

            decomposicao_benezo(N, dNdV, V, T);

            // Euler integration
            for (size_t i = 0; i < N.size(); ++i) {
                N[i] += dV * dNdV[i];
            }

            result.push_back(N);
        }

        // Convert result to MATLAB array
        ArrayFactory factory;
        TypedArray<double> output = factory.createArray<double>({result.size(), 4});
        for (size_t i = 0; i < result.size(); ++i) {
            for (size_t j = 0; j < 4; ++j) {
                output[i][j] = result[i][j];
            }
        }

        outputs[0] = output;
    }

private:
void checkArguments(const ArgumentList& inputs) {
    auto& non_const_inputs = const_cast<ArgumentList&>(inputs);
    if (non_const_inputs.size() != 0) {
        throw std::invalid_argument("This function does not take any input arguments.");
    }
}
};