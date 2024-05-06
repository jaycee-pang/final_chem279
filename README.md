# final_chem279
Chem279 Final Project
This repository contains code for Cholesky decomposition algorithms (with pivoting), testing and comparison of various decomposition algorithms, and applications of Cholesky decomposition such as solving systems of linear equations and constructing inverse matrices. 

To generate the object files and run the executables: 
```bash
make all
```
To run an exectuable: 
```bash
./executable
```
To remove the object files and any txt files not moved into the data directory: 
```bash
make clean
```

To move the generated results files to the data directory: 
```bash
make move
```

To remove the txt files in the data directory: 
```bash
make rmtxt
```

# Directory 
The source files containing the algorithms and testing functions include: 
- `cholesky.cpp`: contains my implementations of Cholesky decomposition of a matrix A with partial pivoting, Cholesky with full pivoting, and LU decomposition. 
- `matgen.cpp`: contains functions for generating matrices that are symmetric positive definite or nearly singular 
- `testing.cpp `: is like the testing 'suite' of functions for applications of cholesky decomposition such as solving systems of linear equations, solving for matrix inverses, recording runtimes, and recording matrix reconstruction errors 

Executables: 
- `test_cholesky`: this was my initial tests for small matrices. Useful for seeing the A, factors L and Lt, and reconstructed A matrices. 
- `time_chol`: run timed tests for the decomposition. Write elapsed time to .txt files for increasing matrix sizes. 
- `chol_err`: run tests to find reconstruction error of original matrix and reconstructed matrix, write results to file for error for increasing matrix size. Write results to txt files 
- `solvers`: teseting applications of Cholesky such as linear solve and inverse matrix problems 
- `test_other`: don't need to run, this is just another executable I used to test some of my implementations


data 
- this directory is where the results for the tests are moved to with the command 
```
make move
```
- `plots.ipynb`: contains a jupyter notebook for reading in the data files and plotting the results 
- plots (sub-directory): contains selected .png files of the plots I included in the report 







