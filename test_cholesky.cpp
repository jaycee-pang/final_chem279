#include <iostream> 
#include "cholesky.h"
#include "matgen.h"
#include <exception> 

/*
Initial testing for small matrices for Cholesky, pivoted Cholesky, LU decomposition, LU with pivoting, 
full pivoted Cholesky, also testing Cholesky with nearly singular matrices. 
*/
int main(void) {
    arma::Mat<double> A = {
        {4,  12, -16},
        {12, 37, -43},
        {-16, -43, 98}
    }; 
    // arma::mat A = {{4,1,-2}, {1,5, 3}, {-2,3,6}};
    std::cout << A.is_sympd() << std::endl;
    try {
        A.print("A");
        std::pair<arma::Mat<double> , arma::Mat<double>> result = pivoted_cholesky(A,false);
        std::cout <<std::endl;
        
        std::cout << "\n L:"<< std::endl;
        result.first.print();  
        std::cout << "\n L trnaspose:"<< std::endl; // upper triangular
        result.second.print(); 
        arma::Mat<double> reconstructedA = result.first*result.second; 
        reconstructedA.print("reconstructed"); 
        A.print("A again");
    }

    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    
    std::cout << "\n\nPivoted Cholesky" << std::endl;

    arma::Mat<double> A2 = gen_sympd(5);
    // arma::mat A2 = {
    //     {4,  12, -16},
    //     {12, 37, -43},
    //     {-16, -43, 98}
    // }; 
    std::cout << A2.is_sympd() << std::endl;
    try {
        A2.print("A");
        std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A2,true);  
        std::cout <<std::endl;
        std::cout << "\n L:"<< std::endl;
        result.first.print(); 
        std::cout << "\n L.t:"<< std::endl;
        result.second.print(); 
        arma::Mat<double> reconstructedA = result.first*result.second; 
        reconstructedA.print("reconstructed"); 
        A2.print("A again");

    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    std::cout << "\n\nArmadillo Cholesky" << std::endl;

    // arma::mat A3 =  {
    //     {4,  12, -16},
    //     {12, 37, -43},
    //     {-16, -43, 98}
    // };
    arma::Mat<double> A3 = gen_sympd(5);
    std::cout << A3.is_sympd() << std::endl;
    try {
        A3.print("A");
        arma::mat L = arma::chol(A3);
        L.print("L"); 
        arma::mat Lt = L.t();
        Lt.print("Lt"); 

        arma::Mat<double> reconstructedA = L * Lt; 
        reconstructedA.print("reconstructed"); 

    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    std::cout << "\n\nLU decomposition" << std::endl;
    arma::Mat<double> A6 = gen_sympd(5);
    try {
        A6.print("A");
        std::pair<arma::Mat<double> , arma::Mat<double>> result = LU_decomp(A6,false);
        std::cout <<std::endl;
        
        std::cout << "\n L:"<< std::endl;
        result.first.print();  
        std::cout << "\n U:"<< std::endl; // upper triangular
        result.second.print(); 
        arma::Mat<double> reconstructedA = result.first*result.second; 
        reconstructedA.print("reconstructed"); 
        
        
    
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }


    std::cout << "\n\nLU decomposition with pivoting" << std::endl;
    arma::Mat<double> A7 = gen_sympd(5);
    try {
        A7.print("A");
        std::pair<arma::Mat<double> , arma::Mat<double>> result = LU_decomp(A7,true);
        std::cout <<std::endl;
        
        std::cout << "\n L:"<< std::endl;
        result.first.print();  
        std::cout << "\n U:"<< std::endl; // upper triangular
        result.second.print(); 
        arma::Mat<double> reconstructedA = result.first*result.second; 
        reconstructedA.print("reconstructed"); 
    
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    // TESTING NEARLY SINGULAR MATRICES 
    arma::Mat<double> A8 = gen_singular(5); 
    std::cout << "\n\n\n Testing nearly singular matrices" <<std::endl;
    try {
        A8.print("A");
        std::pair<arma::Mat<double>, arma::Mat<double>> result = pivoted_cholesky(A8,true);  
        std::cout <<std::endl;
        std::cout << "\n L:"<< std::endl;
        result.first.print(); 
        std::cout << "\n L.t:"<< std::endl;
        result.second.print(); 
        arma::Mat<double> reconstructedA = result.first*result.second; 
        reconstructedA.print("reconstructed"); 
    
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }


    std::cout << "\n\nFull pivoted Cholesky" << std::endl;
    arma::Mat<double> A10 = gen_sympd(5);
    try {
        A10.print("A");
        std::pair<arma::Mat<double>, arma::Mat<double>> result3 = full_pivoted_cholesky(A10);
        std::cout <<std::endl;
        
        std::cout << "\n L:"<< std::endl;
        result3.first.print();  
        std::cout << "\n L trnaspose:"<< std::endl; // upper triangular
        result3.second.print(); 
        arma::Mat<double> reconstructedA = result3.first*result3.second; 
        reconstructedA.print("reconstructed"); 
       

    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    

    arma::Mat<double> A11 = gen_singular(5); 
    std::cout << "\n\n\n Testing nearly singular matrices with full pivoted Chol" <<std::endl;
    try {
        A8.print("A");
        std::pair<arma::Mat<double>, arma::Mat<double>> result4 = full_pivoted_cholesky(A11);  
        std::cout <<std::endl;
        std::cout << "\n L:"<< std::endl;
        result4.first.print(); 
        std::cout << "\n L.t:"<< std::endl;
        result4.second.print(); 
        arma::Mat<double> reconstructedA = result4.first*result4.second; 
        reconstructedA.print("reconstructed"); 
    
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }



    
    return 0; 
}