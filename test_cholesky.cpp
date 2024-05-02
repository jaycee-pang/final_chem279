#include <iostream> 
#include "cholesky.h"
#include "matgen.h"
#include <exception> 


int main(void) {
    // arma::mat A = gen_sympd(5);
    arma::mat A = {
        {4,  12, -16},
        {12, 37, -43},
        {-16, -43, 98}
    }; 
    // arma::mat A = {{4,1,-2}, {1,5, 3}, {-2,3,6}};
    
    std::cout << A.is_sympd() << std::endl;

    try {
        A.print("A");
        std::pair<arma::mat , arma::mat> result = cholesky(A);
        std::cout <<std::endl;
        
        std::cout << "\n L:"<< std::endl;
        result.first.print();  
        std::cout << "\n L trnaspose:"<< std::endl; // upper triangular
        result.second.print(); 
        arma::mat reconstructedA = result.first*result.second; 
        reconstructedA.print("reconstructed"); 
       

    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    
    std::cout << "\n\nPivoted Cholesky" << std::endl;

    arma::mat A2 = {{4,1,-2}, {1,5, 3}, {-2,3,6}};
    std::cout << A2.is_sympd() << std::endl;
    try {
        A2.print("A");
        std::pair<arma::mat, arma::mat> result = pivoted_cholesky(A,true);  
        std::cout <<std::endl;
        std::cout << "\n L:"<< std::endl;
        result.first.print(); 
        std::cout << "\n L.t:"<< std::endl;
        result.second.print(); 
        arma::mat reconstructedA = result.first*result.second; 
        reconstructedA.print("reconstructed"); 

    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    std::cout << "\n\nArmadillo Cholesky" << std::endl;

    arma::mat A3 =  {
        {4,  12, -16},
        {12, 37, -43},
        {-16, -43, 98}
    };
    std::cout << A3.is_sympd() << std::endl;
    try {
        A3.print("A");
        arma::mat L = arma::chol(A);
        L.print("L"); 
        arma::mat Lt = L.t();
        Lt.print("Lt"); 

        arma::mat reconstructedA = L * Lt; 
        reconstructedA.print("reconstructed"); 

    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }




    // std::cout << "\n\nNew trial chol" << std::endl;
    // A = {
    //     {4,  12, -16},
    //     {12, 37, -43},
    //     {-16, -43, 98}
    // };
    // try {
    //     A.print("A");
    //     std::pair<arma::mat , arma::mat> result = other_chol(A,true);
    //     std::cout <<std::endl;
        
    //     std::cout << "\n L:"<< std::endl;
    //     result.first.print();  
    //     std::cout << "\n L trnaspose:"<< std::endl; // upper triangular
    //     result.second.print(); 
    //     arma::mat reconstructedA = result.first*result.second; 
    //     reconstructedA.print("reconstructed"); 
       

    // }
    // catch (const std::exception& e) {
    //     std::cerr << e.what() << std::endl;
    // }





    
    return 0; 
}