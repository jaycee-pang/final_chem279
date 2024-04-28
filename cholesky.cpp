#include "cholesky.cpp"
/*
.is_trimatu / .is_trimatl	 	check whether matrix is upper/lower triangular
.is_diagmat	 	check whether matrix is diagonal
.is_square	 	check whether matrix is square sized
.is_symmetric	 	check whether matrix is symmetric
.is_hermitian	 	check whether matrix is hermitian
.is_sympd	 	check whether matrix is symmetric/hermitian positive definite
*/

int find_pivot(const mat& A, int start) {
    int n = A.n_rows();
    int pivot_row = start; 
    double max_val = std::abs(A(start, start)); 
    // find the max val 
    for (int i=start+1; i<n ; i++) {
        if (std::abs(A(i,start)) > max_val) {
            max_val = std::abs(A(i, start));
            pivot_row = i; 

        }
        
    }
    return pivot_row;
}

/*
A  = LL.t 
diagonals: Lkk= sqrt(Akk- sum from j=1 to k-1 of [Lkj^2])
off-diagonals for i>k : Lik = (1/Lkk)(Aik - sum from j=1 to k-1 of [Lij*Lkj])
return lower cholesky factor 
pivoted: 
(P.t)AP=L(L.t)

L is lower triangular
*/
bool cholesky(const arma::mat & A, arma::mat& L, bool pivot) {
    int n = A.n_rows();
    L.zeros(); 
    arma::uvec p(n); 
    for (int i=0; i<n; i++) {
        for (int j=0; j<=i; j++) {
            double sum=0.0; 
            for (int k=0; k<j; k++) {
                sum += L(i, k)*L(j,k);

                if (i==j) { 
                    // check for positive definiteness 
                    if (A(i,i) - sum <= 0) {return false; }
                    L(i, j) = std::sqrt(A(i,i) - sum);
                }
                else {
                    if (L(j,j) == 0) {return false;}
                    L(i, j) = (1.0 /L(j, j))*(A(i, j) - sum);
                }
            }
        }
    }
    return true;

}


bool pivoted_cholesky(const mat& A, mat& L) {
    int n = A.n_rows();
    L.arma::zeros(); 
    arma::uvec p(n);
    p.arma::zeros(); 

    for (int i=0; i<n; i++) {
        int pivot_row = find_pivot(A,i); 
        if (pivot_row != i) {
            // swap rows i, pivot row of A 
            // swap i, pivot rows of L 
            p(i) = pivot_row;

        }

        for (int j=i; j<n; j++) {
            double sum=0.0; 
            for (int k=0; k<i; k++) {
                sum += L(j,k)*L(i, k);

                if (i==j) {
                    // check for positive definiteness 
                    if (A(i,i) - sum <= 0) {return false; }
                    L(i,i) = std::sqrt(A(i,i) - sum);
                }
                else {
                    if (L(j,j) == 0) {return false;} // zero division
                    L(j,i) = (1.0 / L(j,j)) * (A(j,i) - sum);

                }
            
                
                // double diagonal = L(j,j); 
                // for (int k=i+1; k<n; k++) {
                //     L(k,i) /= digonal;
                // }
                
                // L(i, j) = (1.0 /L(j, j))*(A(i, j) - sum);
            
            }
        }

    return true;
    }


        
        
}