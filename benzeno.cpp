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
const double R = 8314.0;
const double P = 101325.0;


typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

std::ofstream myfile("output.csv");

double calculate_k(double A, double E, double T) {
    return (A * (std::exp((-1*E) / (R * T))));
}

double calculate_concentration(double N_species, double N_total, double T) {
    return ((P / (R * T)) * (N_species / N_total));
}

void decomposicao_benezo(const vector_type& N, vector_type& dNdV, double V, double T) {
    double _k1 = 7.0e5;  // L/mol.h
    double K1 = 0.31;
    double _k2 = 4e5; // L/mol.h
    double K2 = 0.48;

    // Fluxos molares
    double N_total = std::accumulate(N.begin(), N.end(), 0.0);

    // Concentrações
    double C_c6h6 = calculate_concentration(N[0], N_total, T); //m3/mol para L/mol
    double C_h = calculate_concentration(N[1], N_total, T);
    double C_c12h10 = calculate_concentration(N[2], N_total, T);
    double C_c18h14 = calculate_concentration(N[3], N_total, T);

    double r1 = _k1 * ( std::pow(C_c6h6, 2) - ((C_c12h10 * C_h)/K1));
    double r2 = _k2 * ((C_c6h6 * C_c12h10)  - ((C_c18h14 * C_h)/K2));

    dNdV[0] =  -(2*r1) - r2;    // dNc2h6/dV
    dNdV[1] =  r1 + r2;         // dNno/dV
    dNdV[2] =  r1 - r2;         // dNc2h5/dV
    dNdV[3] =  r2;              // dNhno/dV
}

struct sysBenzeno {
    double T = 1033; // Temperatura do sistema

    // Método para calcular as derivadas (EDOs)
    void operator()(const vector_type& N, vector_type& dNdV, double V) const {
        decomposicao_benezo(N, dNdV, V, T);
    }
};

struct JBenzeno {
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

        typedef runge_kutta_dopri5< vector_type > dopri5_type;
        // typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
        // typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;

        //Benzeno
        vector_type B0 (4);
        for (int i=0; i < 4; i++) { B0[i] = 0.0; }
        B0[0] = 1.0;

        // size_t num_of_steps = integrate_adaptive(make_controlled(0.1, 0.1,  dopri5_type() )  ,
        //     sysBenzeno(), B0, 0.0, 1500.0, (1500.0/1e9),
        //     myObserver
        // );

        size_t num_of_steps = integrate_const( make_dense_output< rosenbrock4< double > > (1.0e-2, 1.0e-2) ,
            make_pair(sysBenzeno() , JBenzeno()),
            B0 , 0.0 , 1500.0, (1500.0/1.0e6), 
            myObserver
        );

        clog << num_of_steps << endl;

    } catch (std::exception &e) {
        std::cerr << "Erro: " << e.what() <<std:: endl;
    }


    myfile.close();
    return 0;
}
