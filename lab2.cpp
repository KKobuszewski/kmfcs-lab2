#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <complex>
#include <complex.h>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <gaussians.h>


void get_overlap_matrix(const unsigned N, const double alpha, Eigen::MatrixXcd & S_mat)
{
    const double a0 = 1.;
    // iterate over quanatum numbers
    for (unsigned i=0; i<N; i++)
    for (unsigned j=0; j<N; j++)
    for (unsigned k=0; k<N; k++)
    for (unsigned l=0; l<N; l++)
    {
#ifdef DEBUG
        unsigned i_idx = i*N+j;
        unsigned j_idx = k*N+l;
        std::cout << i << " " << j << "\t" << "(" << i_idx << "," << j_idx << ")" << get_S(i+1,j+1,k+1,l+1,alpha,a0) << std::endl;
#endif
        S_mat(i*N+j,k*N+l) = (std::complex<double>) (get_S(i+1,j+1,k+1,l+1,alpha,a0) + I*0.);
    }
}

void get_hamiltonian_matrix(const unsigned N, const double alpha, Eigen::MatrixXcd & H_mat)
{
    const double a0 = 1.;
    // iterate over quanatum numbers
    for (unsigned i=0; i<N; i++)
    for (unsigned j=0; j<N; j++)
    for (unsigned k=0; k<N; k++)
    for (unsigned l=0; l<N; l++)
    {
#ifdef DEBUG
        unsigned i_idx = i*N+j;
        unsigned j_idx = k*N+l;
        std::cout << i << " " << j << "\t" << "(" << i_idx << "," << j_idx << ")" << std::endl;
#endif
        H_mat(i*N+j,k*N+l) = (std::complex<double>) (get_T(i+1,j+1,k+1,l+1,alpha,a0) + get_V(i+1,j+1,k+1,l+1,1,alpha,a0) + get_V(i+1,j+1,k+1,l+1,2,alpha,a0) + I*0.);
    }
}





int main()
{
    double _alpha = 1.0; // 9.5 to 10.5
    unsigned N = 4; // 5 to 10
    
    for (N; N < 10; N++)
    {
        std::cout << "N: " << N << std::endl;
        
        Eigen::MatrixXcd S_mat(N*N,N*N);
        Eigen::MatrixXcd H_mat(N*N,N*N);
        
        std::ostringstream oss;
        oss << "energies_N" << N << ".dat";
        
        std::ofstream file;
        file.open(oss.str().c_str());
        
        if (N>6) _alpha = 5.0;
        if (N>7) _alpha = 35.0;
        if (N>8) _alpha = 55.0;
        if (N>9) _alpha = 75.0;
        
        for (double alpha = _alpha; alpha <= 100.0; alpha += 0.1)
        {
            get_overlap_matrix(N, alpha, S_mat);
#ifdef DEBUG
            std::cout << S_mat << std::endl << std::endl;
#endif
            
            get_hamiltonian_matrix(N, alpha, H_mat);
#ifdef DEBUG
            std::cout << H_mat << std::endl << std::endl;
#endif
            
            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> es(H_mat,S_mat,Eigen::Ax_lBx | Eigen::EigenvaluesOnly);
        
            //std::cout << " Eigenvalues: " << std::endl;
            file << std::setprecision(15) << std::fixed;
            file << alpha << "\t" << es.eigenvalues().transpose() << std::endl;
            std::cout << std::setprecision(9) << std::fixed;
            std::cout << alpha << "\t" << es.eigenvalues()(0) << std::endl;
        }
        std::cout << std::endl;
    }
    
    return EXIT_SUCCESS;
}