#include <iostream>
#include <fstream>
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
        std::vector<double> V_span = {0.0, 1500.0};     // Volume span
        double T = 1033.0;                            // Temperature (K)

        // Create MATLAB data arrays
        TypedArray<double> N0_matlab = matlabPtr->createArray({4}, N0.begin(), N0.end());
        TypedArray<double> V_span_matlab = matlabPtr->createArray({2}, V_span.begin(), V_span.end());

        // Call MATLAB function
        auto results = matlabPtr->feval(u"solveBenzenoODE", {N0_matlab, V_span_matlab, matlab::data::ArrayFactory().createScalar(T)});

        // Extract results
        auto V_result = results[0].getArray<double>();
        auto N_result = results[1].getArray<double>();

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