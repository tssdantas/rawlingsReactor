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
const double A1 = 1.0e14;
const double E1 = 217600.0;
const double A2 = 3.0e14; 
const double E2 = 165300.0;
const double A3 = 3.4e12; 
const double E3 = 28500.0;
const double A4 = 1.0e12; 
const double E4 = 0.0;
const double A5 = 1.0e13; 
const double E5 = 200800.0;
const double A6 = 1.0e12; 
const double E6 = 0.0;
const double R = 8.314;
const double P = 101325.0;


typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

double calculate_k(double A, double E, double T) {
    return A * std::exp(-E / (R * T));
}

double calculate_concentration(double N_species, double N_total, double T) {
    return (P / (R * T)) * (N_species / N_total);
}

// double calculate_concentration_P(double N_species, double N_total, double T) {
//     return ((R * T) / P) * (N_species / N_total);
// }

// Função do sistema de equações diferenciais
void decomposicao_etano(const vector_type& N, vector_type& dNdV, double V, double T) {
    // Calculando as constantes de taxa
    double k1 = calculate_k(A1, E1, T);
    double k2 = calculate_k(A2, E2, T);
    double k3 = calculate_k(A3, E3, T);
    double k4 = calculate_k(A4, E4, T);
    double k5 = calculate_k(A5, E5, T);
    double k6 = calculate_k(A6, E6, T);

    // Fluxos molares
    double N_total = std::accumulate(N.begin(), N.end(), 0.0);

    // Concentrações
    double C_c2h6 = calculate_concentration(N[0], N_total, T);
    double C_no = calculate_concentration(N[1], N_total, T);
    double C_c2h5 = calculate_concentration(N[2], N_total, T);
    double C_hno = calculate_concentration(N[3], N_total, T);
    double C_h = calculate_concentration(N[4], N_total, T);

    // Taxas de reação
    double r1 = k1 * C_c2h6 * C_no;
    double r2 = k2 * C_c2h5;
    double r3 = k3 * C_h * C_c2h6;
    double r4 = k4 * C_h * C_no;
    double r5 = k5 * C_hno;
    double r6 = k6 * C_c2h5 * C_hno;

    // Equações diferenciais
    dNdV[0] =  -r1 - r3 + r6;        // dNc2h6/dV
    dNdV[1] =  -r1 - r4 + r5 + r6;   // dNno/dV
    dNdV[2] =  r1 - r2 + r3 - r6;    // dNc2h5/dV
    dNdV[3] =  r1 + r4 - r5 - r6;    // dNhno/dV
    dNdV[4] =  r2 - r3 - r4 + r5;    // dNh/dV
    dNdV[5] =  r2;                   // dNc2h4/dV
    dNdV[6] =  r3;                   // dNh2/dV
}

void decomposicao_benezo(const vector_type& N, vector_type& dNdV, double V, double T) {
    double _k1 = 7.0e5;
    double K1 = 0.31;
    double _k2 = 4e5;
    double K2 = 0.48;

    // Fluxos molares
    double N_total = std::accumulate(N.begin(), N.end(), 0.0);

    // Concentrações
    double C_c6h6 = calculate_concentration(N[0], N_total, T);
    double C_h = calculate_concentration(N[1], N_total, T);
    double C_c12h10 = calculate_concentration(N[2], N_total, T);
    double C_c18h14 = calculate_concentration(N[3], N_total, T);

    double r1 = _k1 * ( std::pow(C_c6h6, 2) - ((C_c12h10 * C_h)/K1));
    double r2 = _k2 * ((C_c6h6 * C_c12h10)  - ((C_c18h14 * C_h)/K2));

    dNdV[0] =  -(2*r1) - r2;    // dNc2h6/dV
    dNdV[1] =  r1 + r2;         // dNno/dV
    dNdV[2] =  r1 - r2;    // dNc2h5/dV
    dNdV[3] =  r2;    // dNhno/dV
}

struct sysBenzeno {
    double T = 1000; // Temperatura do sistema

    // Método para calcular as derivadas (EDOs)
    void operator()(const vector_type& N, vector_type& dNdV, double V) const {
        decomposicao_benezo(N, dNdV, V, T);
    }
};

struct sysEtano {
    double T = 1000; // Temperatura do sistema

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
        for (auto _v:x) {
            std::cout << _v << "\t";
        }
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
        N0[1] = 3.48e-4;

        typedef runge_kutta_dopri5< vector_type > dopri5_type;
        // typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
        // typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;

        // size_t num_of_steps = integrate_const( make_dense_output( 3.0e-3 , 2.0e-3 , dopri5_type() ),
        //     sysEtano(), N0, 0.0 , 0.0015, (0.0015/1e3),
        //     myObserver
        // );


        // size_t num_of_steps = integrate_n_steps( make_dense_output( 1.0e-2 , 1.0e-2 ,  runge_kutta_dopri5< vector_type >()) ,
        //     sysEtano(), N0, 0.0 , 0.0015 , (0.0015/1e2),
        //     myObserver
        // );


        // size_t num_of_steps = integrate_adaptive(make_controlled( 1.0e-3 , 1.0e-3 ,  dopri5_type() )  ,
        //     sysEtano(), N0, 0.0, 0.0015, (0.0015/1e6),
        //     myObserver
        // );

        // size_t num_of_steps = integrate_const( make_dense_output< rosenbrock4< double > > (5.0e-3, 5.0e-3) ,
        //         make_pair(sysEtano() , JEtano()) ,
        //         N0 , 0.0 , 0.0015, (0.0015/1e3), 
        //         myObserver
        // );
        
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

        //Benzeno
        vector_type B0 (4);
        for (int i=0; i < 4; i++) { B0[i] = 0.0; }
        B0[0] = 1.0;

        size_t num_of_steps = integrate_adaptive(make_controlled( 1.0e-3 , 1.0e-3 ,  dopri5_type() )  ,
            sysBenzeno(), B0, 0.0, 0.0015, (0.0015/1e6),
            myObserver
        );

        clog << num_of_steps << endl;

    } catch (std::exception &e) {
        std::cerr << "Erro: " << e.what() <<std:: endl;
    }



    return 0;
}
