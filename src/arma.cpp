#include <vector> 
#include <algorithm>
#include <stdexcept>  

// matrix should be real, 2D only
class mat {
    private: 
        int m_; // number of rows 
        int n_; // number of cols 
        std::vector<std::vector<double>> mat_; // vector of vectors for now arr is faster
    public: 
        // default constructor 
        mat(): m_(0), n_(0), mat(0, std::vector<double>(0)) {
            mat_.resize(m_, std::vector<double>(0)); 
        }
        // initialize size 
        mat(int m, int n): m_(m), n_(n), mat_(m, std::vector<double>(n)) {
            mat_.resize(m, std::vector<double>(n));
        }
        // initialize from another matrix (vector of vectors )
        mat(std::vector<std::vector<double>> & matrix): m_(matrix.size()), n_(matrix[0].size(), mat_(matrix)) {} 
        // copy constructor 
        mat(const mat& other): m_(other.m_), n_(other.n_), mat_(other.mat_) {}

        // getters 
        int m() const {return m_;}
        int n() const {return n_;}
        const std::vector<std::vector<double>>& mat() const {return mat_;}

        void reshape(int m, int n) {
            m_ = m; 
            n_= n;
            mat_.resize(m, std::vector<double>(n));
        }
        void print() const {
            for (const auto &m : mat_) {
                for (const auto& n: m) {
                    std::cout << n << " "; 
                }
                std::cout << std::endl; 
            }
        }


        // operators 
        // non-const to modify via = 
        const double operator()(int row, int col) {
            return mat_[row][col]; 
        }
        // const ver, return reference 
        const double& operator()(int row, int col) const {
            return mat_[row][col]; 
        }
        /*
        operators: +, -, ==, *, 
        scalar operations 
        transpose (in place? )

        */



}; 

class 