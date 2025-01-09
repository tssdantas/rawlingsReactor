#include <iostream>
#include <vector>
#include <numeric>
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

using namespace matlab::engine;
using namespace matlab::data;

int main(int argc, char** argv) {
    try {
        // Start MATLAB session
        std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();

        // Define initial conditions and parameters
        std::vector<double> N0 = {1.0, 0.0, 0.0, 0.0}; // Initial moles
        std::vector<double> V_span = {0.0, 10.0};     // Volume span
        double T = 1033.0;                            // Temperature (K)

        // Use ArrayFactory to create MATLAB arrays
        ArrayFactory factory;
        TypedArray<double> N0_matlab = factory.createArray<double>({4}, N0.data(), N0.data() + N0.size());
        TypedArray<double> V_span_matlab = factory.createArray<double>({2}, V_span.data(), V_span.data() + V_span.size());

        // Call MATLAB function
        std::vector<Array> args = {N0_matlab, V_span_matlab, factory.createScalar(T)};
        std::vector<Array> results = matlabPtr->feval(u"solveBenzenoODE", 2, args);

        // Extract results
        TypedArray<double> V_result = results[0];
        TypedArray<double> N_result = results[1];

        // Print the results
        std::cout << "Volume (V): ";
        for (const auto& V : V_result) {
            std::cout << V << " ";
        }
        std::cout << "\n";

        std::cout << "Mole Fractions (N): ";
        for (const auto& N : N_result) {
            std::cout << N << " ";
        }
        std::cout << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}