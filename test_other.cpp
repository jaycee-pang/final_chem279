#include "cholesky.h"
#include "matgen.h"
#include "testing.h"
#include <iostream> 
// testing 

int main(void) {
    std::vector<int> mat_sizes = {3,4,5}; 
    std::vector<arma::Mat<double>> matrices;
    // so we can test the same matrices for each method/algo.
    for (int size : mat_sizes) {
        matrices.push_back(gen_sympd(size));
    }
    std::cout << "Cholesky matrix reconstruction error" << std::endl;
    
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i]; 
        double error1 = chol_err(A, false); 
    
        
        
    }

    std::cout << "\n\nPivoted Cholesky matrix reconstruction error" << std::endl;

    
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A2 = matrices[i];
        double error2 = chol_err(A2, true); 
 
        
    }


  
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A = matrices[i];
        arma::Mat<double> L = arma::chol(A); 
        arma::Mat<double> reconstructedA3 = L*L.t(); 
        arma::Mat<double> diff3 = A - reconstructedA3; 
        double error3 = arma::norm(diff3, "fro"); 
  
        
    }



    std::cout << "\n\nPivoted Cholesky matrix reconstruction error on nearly singular matrices" << std::endl;
    
    std::vector<arma::Mat<double>> almost_singular;
    for (int size : mat_sizes) {
        almost_singular.push_back(gen_singular(size));
    }
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A4 = almost_singular[i]; 
        double error4 = chol_err(A4, true); 

        
        
    }

    std::cout << "\n\nFull Pivoted Cholesky matrix reconstruction error" << std::endl;
    
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A5 = matrices[i]; 
        // A5.print("A");
        arma::Mat<double> og_A5 = A5; 
        std::pair<arma::Mat<double>, arma::Mat<double>> result5 = full_pivoted_cholesky(A5); 
        // double error = chol_err(A, false); 
        arma::Mat<double> L5= result5.first; 
        arma::Mat<double> Lt5 = result5.second; 
        arma::Mat<double> reconA5 = L5*Lt5;
        // reconA5.print("reconstructed A"); 
        arma::Mat<double> diff5 =og_A5 -  reconA5; 
        double error5 = arma::norm(diff5, "fro"); 
    

        
    }
 

    std::cout << "\n\nPivoted Cholesky matrix reconstruction error on nearly singular matrices" << std::endl;
   
    for (int i=0; i<mat_sizes.size(); i++) {
        int size = mat_sizes[i];
        arma::Mat<double> A6 = almost_singular[i]; 
        // A6.print("A"); 
        arma::Mat<double> og_A6 = A6; 
        std::pair<arma::Mat<double>, arma::Mat<double>> result6 = full_pivoted_cholesky(A6); 
        arma::Mat<double> L6 = result6.first; 
        arma::Mat<double> Lt6 = result6.second; 
        arma::Mat<double> reconA6 = L6*Lt6;
        arma::Mat<double> diff6 =og_A6 -  reconA6; 
        // reconA6.print("reconstructed");
        double error6 = arma::norm(diff6, "fro"); 
 
        
    }


    return 0; 

}
    

