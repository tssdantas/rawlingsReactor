#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>  
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

using vector_type = std::vector<double>;

const double RG2 = 8.314;
const double P = 1.0;

double calculate_k(double A, double E, double T) {
    return A * exp(-E / (RG2 * T));
}

double calculate_concentration(double N, double Q) {
    return N / Q;
}

void decomposicao_etano(const vector_type& N, vector_type& dNdV, double V, double T) {
    const double A11 = 1e14;
    const double A12 = 1e12;
    const double A2 = 3e14; 
    const double A3 = 3.4e12; 
    const double A41 = 1e12; 
    const double A42 = 1e13;

    const double E11 = 217.6e3; 
    const double E12 = 0; 
    const double E2 = 165.3e3; 
    const double E3 = 28.5e3; 
    const double E41 = 0; 
    const double E42 = 200.8e3;

    double k11 = calculate_k(A11, E11, T);
    double k12 = calculate_k(A12, E12, T);
    double k2 = calculate_k(A2, E2, T);
    double k3 = calculate_k(A3, E3, T);
    double k41 = calculate_k(A41, E41, T);
    double k42 = calculate_k(A42, E42, T);

    double N_total = std::accumulate(N.begin(), N.end(), 0.0);
    double Q = ((RG2 * T) / P) * (N_total);

    double C_c2h6 = calculate_concentration(N[0], Q);
    double C_c2h5_p = calculate_concentration(N[1], Q);
    double C_c2h4 = calculate_concentration(N[2], Q);
    double C_h = calculate_concentration(N[3], Q);
    double C_h2 = calculate_concentration(N[4], Q);
    double C_NO = calculate_concentration(N[5], Q);
    double C_HNO = calculate_concentration(N[6], Q);

    double r1 = (k11 * C_c2h6 * C_NO) - (k12 * C_c2h5_p * C_HNO);
    double r2 = k2 * C_c2h5_p;
    double r3 = k3 * C_h * C_c2h6;
    double r4 = (k41 * C_h * C_NO) - (k42 * C_HNO);

    dNdV[0] = (-1) * (r1 + r3);
    dNdV[1] = r1 - r2 + r3;
    dNdV[2] = r2;
    dNdV[3] = r2 - r3 - r4;
    dNdV[4] = r3;
    dNdV[5] = (-1) * (r1 + r4);
    dNdV[6] = r1 + r4;
}

int f(realtype V, N_Vector y, N_Vector ydot, void *user_data) {
    double T = 1050.0; 
    
    vector_type N(7); 
    vector_type dNdV(7);
    for (int i = 0; i < 7; i++) {
        N[i] = NV_Ith_S(y, i);
    }

    decomposicao_etano(N, dNdV, V, T);

    for (int i = 0; i < 7; i++) {
        NV_Ith_S(ydot, i) = dNdV[i];
    }

    return 0;
}

int main() {
    realtype V0 = 0.0;
    realtype V1 = 1500.0;
    realtype reltol = 1e-4;
    realtype abstol = 1e-8;

    N_Vector y = N_VNew_Serial(7);
    NV_Ith_S(y, 0) = 6.62e-3;
    NV_Ith_S(y, 5) = 3.48e-4;
    
    for (int i = 1; i < 7; i++) {
        if (i != 5) NV_Ith_S(y, i) = 0.0;
    }

    void *cvode_mem = CVodeCreate(CV_BDF);
    CVodeInit(cvode_mem, f, V0, y);
    CVodeSStolerances(cvode_mem, reltol, abstol);
    CVodeSetUserData(cvode_mem, NULL);

    SUNMatrix A = SUNDenseMatrix(7, 7);
    SUNLinearSolver LS = SUNLinSol_Dense(y, A);
    CVodeSetLinearSolver(cvode_mem, LS, A);

    std::ofstream file("resultado.csv");
    //file << "V,N0,N1,N2,N3,N4,N5,N6\n";

    realtype V = V0;
    realtype Vstep = 10.0; 

    while (V < V1) {
        CVode(cvode_mem, V + Vstep, y, &V, CV_NORMAL);

        std::cout << "V = " << V;
        file << V; 

        for (int j = 0; j < 7; j++) {
            std::cout << ", N[" << j << "] = " << NV_Ith_S(y, j);
            file << "," << NV_Ith_S(y, j);
        }

        std::cout << std::endl;
        file << "\n";
    }

    file.close();

    N_VDestroy(y);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);

    return 0;
}
