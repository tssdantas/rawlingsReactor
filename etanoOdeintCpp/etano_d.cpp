#include <iostream>
#include <fstream>
#include <utility>
#include <cmath>
#include <boost/numeric/odeint.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;

const double mu = 1000.0;

// Definições de constantes
// const double A1 = 1e14;
// const double E1 = 217600.0;
// const double A2 = 3e14; 
// const double E2 = 165300.0;
// const double A3 = 3.4e12; 
// const double E3 = 28500.0;
// const double A4 = 1e12; 
// const double E4 = 0.0;
// const double A5 = 1e13; 
// const double E5 = 200800.0;
// const double A6 = 1e12; 
// const double E6 = 0.0;

const double A11 = 1e14;
const double A12 = 1e12;
const double A2 = 3e14; 
const double A3=3.4e12; 
const double A41=1e12; 
const double A42=1e13;

const double E11 = 217.6e3; 
const double E12=0; 
const double E2=165.3e3; 
const double E3=28.5e3; 
const double E41=0; 
const double E42=200.8e3;

//Constantes Físicas
const double RG1 = 8.314; //j/mol K
const double RG2 = 82.06; //cm³ atm/ mol K

//Constantes de operação
const double T = 1050;  //K
const double P = 1;     //atm 
const double Qf = 600;  //cm^3/s




typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

std::ofstream myfile("output.csv");

double calculate_k(double A, double E, double T) {
    return (A * (std::exp((-1*E) / (RG1 * T))));
}

double calculate_concentration(double N_species, double _Q) {
    return (N_species / _Q);
}


// Função do sistema de equações diferenciais
void decomposicao_etano(const vector_type& N, vector_type& dNdV, double V, double T) {
    // Calculando as constantes de taxa
    double k11 = calculate_k(A11, E11, T);
    double k12 = calculate_k(A12, E12, T);
    double k2 = calculate_k(A2, E2, T);
    double k3 = calculate_k(A3, E3, T);
    double k41 = calculate_k(A41, E41, T);
    double k42 = calculate_k(A42, E42, T);

    // Fluxos molares
    double N_total = std::accumulate(N.begin(), N.end(), 0.0);

    double Q = ((RG2*T)/P)*(N_total);

    // Concentrações
    double C_c2h6 = calculate_concentration(N[0], Q); //1
    double C_c2h5_p = calculate_concentration(N[1], Q); //2
    double C_c2h4 = calculate_concentration(N[2], Q); //3
    double C_h = calculate_concentration(N[3], Q); //4
    double C_h2 = calculate_concentration(N[4], Q); //5 
    double C_NO = calculate_concentration(N[5], Q); //6
    double C_HNO = calculate_concentration(N[6], Q); //7

    // Taxas de reação
    double r1 = (k11 * C_c2h6 * C_NO) - (k12 * C_c2h5_p * C_HNO);
    double r2 = k2 * C_c2h5_p;
    double r3 = k3 * C_h * C_c2h6;
    double r4 = (k41 * C_h * C_NO) - (k42*C_HNO);
    // double r5 = k5 * C_hno;
    // double r6 = k6 * C_c2h5 * C_hno;

    // Equações diferenciais
    dNdV[0] =  (-1) * (r1 + r3);        // dNc2h6/dV
    dNdV[1] =  r1 - r2 + r3;   // dNno/dV
    dNdV[2] =  r2;    // dNc2h5/dV
    dNdV[3] =  r2 - r3 - r4;    // dNhno/dV
    dNdV[4] =  r3;    // dNh/dV
    dNdV[5] =  (-1) * (r1 + r4) ;                   // dNc2h4/dV
    dNdV[6] =  r1+r4;                   // dNh2/dV
}

struct sysEtano {
    double T = 1050; // Temperatura do sistema

    // Método para calcular as derivadas (EDOs)
    void operator()(const vector_type& N, vector_type& dNdV, double V) const {
        decomposicao_etano(N, dNdV, V, T);
    }
};

struct JEtano {
    double T = 1050; // Temperatura do sistema

    // Método para calcular as derivadas (EDOs)
    void operator()(const vector_type &N , matrix_type &J , const double &V , vector_type &dNdV) const {
        //
        J.clear();
    }
};

void myObserver(const vector_type &x, double t) {
    try {
        std::cout << t << '\t';
        myfile << t << ", "; 
        for (auto _v:x) {
            myfile << _v << ", ";
            std::cout << _v << "\t";
        }
        myfile << std::endl;
        std::cout << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "Erro: " << e.what() <<std:: endl;
    }

}

struct StepOverflowChecker {
    size_t max_steps;
    size_t current_step;

    StepOverflowChecker(size_t max_steps_) : max_steps(max_steps_), current_step(0) {}

    template <class State>
    void operator()(const State &, double) {
        if (++current_step > max_steps) {
            throw std::overflow_error("Maximum number of steps exceeded!");
        }
    }
};



int main( int argc , char **argv )
{
    /* initialize random seed: */
    srand ( time(NULL) );

    try {

        // Condições iniciais Etano
        vector_type N0 (7);
        for (int i=0; i<7; i++) { N0[i] = 0.0; }

        N0[0] = 6.62e-3;
        N0[5] = 3.48e-4;

        N0[0] = ((0.95*Qf*P)/(RG2*T));
        N0[5] = ((0.05*Qf*P)/(RG2*T));

        typedef runge_kutta_dopri5< vector_type > dopri5_type;
        // typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
        // typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;

        // size_t num_of_steps = integrate_const( make_dense_output( 1.0e-2 , 1.0e-2 , dopri5_type() ),
        //     sysEtano(), N0, 0.0 , 1500.0, (1500.0/1e3),
        //     myObserver
        // );


        // size_t num_of_steps = integrate_n_steps( make_dense_output( 1.0e-1, 1.0e-1, 1.0e-1,  dopri5_type()) ,
        //     sysEtano(), N0, 0.0 , (0.0015/1e3), 100,
        //     myObserver
        // );


        // size_t num_of_steps = integrate_adaptive(make_controlled( 1.0e-2 , 1.0e-2 ,  dopri5_type() )  ,
        //     sysEtano(), N0, 0.0, 1500.0, (1500/1e4),
        //     myObserver
        // );

        size_t num_of_steps = integrate_const( make_dense_output< rosenbrock4< double > > (1.0e-1, 1.0e-1) ,
                make_pair(sysEtano() , JEtano()) ,
                N0 , 0.0 , 1500.0, (1500/1e3), 
                myObserver
        );
        
        //clog << num_of_steps << endl;

        //int steps = 1000;
        // StepOverflowChecker overflow_checker(steps++); // Maximum 100 steps allowed

        // size_t num_of_steps_A = integrate_n_steps(
        //     make_dense_output( 1.0e-2 , 1.0e-2 ,  runge_kutta_dopri5< vector_type >()),
        //     sysEtano(),    // System
        //     N0,            // Initial state
        //     0.0,            // Start time
        //     (1500.0 / 1e4), // Time step
        //     steps,           // Number of steps
        //     [&overflow_checker](const vector_type &N, double t) {
        //         overflow_checker(N, t); // Check for step overflow
        //         myObserver(N, t);     // Call observer
        //     }
        // );


        clog << num_of_steps << endl;

    } catch (std::exception &e) {
        std::cerr << "Erro: " << e.what() <<std:: endl;
    }


    myfile.close();
    return 0;
}
