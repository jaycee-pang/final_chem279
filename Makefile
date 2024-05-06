CPP = g++
C = gcc
CPPFLAGS = -std=c++20 -I/opt/homebrew/opt/armadillo/include -O1 -w
INCLUDES = -Iinc -I/opt/homebrew/opt/armadillo/include
LDFLAGS = -L/opt/homebrew/opt/armadillo/lib -larmadillo -lblas -llapack
OBJS = cholesky.o matgen.o testing.o
TXT = pivoted_cholesky_times.txt cholesky_times.txt armachol_times.txt cholesky_error.txt pivoted_cholesky_error.txt \
	arma_error.txt LU_times.txt LU_pivot_times.txt cholesky_inv_err.txt piv_cholesky_inv_err.txt \
	LU_decomp_error.txt LU_pivot_decomp_error.txt pivoted_cholesky_times_sing.txt pivoted_cholesky_err_sing.txt \
	cholesky_Axb_err.txt cholesky_piv_Axb_err.txt full_piv_cholesky_error.txt full_pivoted_cholesky_times.txt \
	cholesky_full_piv_Axb_err.txt full_piv_cholesky_inv_err.txt full_pivoted_cholesky_err_sing.txt full_pivoted_cholesky_sing_times.txt

EXECU = test_cholesky time_chol chol_err solvers test_other

# contains cholesky algorithms 
cholesky.o: cholesky.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c cholesky.cpp -o cholesky.o

# code for generating symmetric positive definite matrices and nearly singular matrices 
matgen.o: matgen.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c matgen.cpp -o matgen.o

# testing functions: run Cholesky and record times, error, etc. 
testing.o: testing.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c testing.cpp -o testing.o

# first tests, tests small matrices for different algorithms 
test_cholesky: $(OBJS)
	$(CPP) $(CPPFLAGS) -o test_cholesky test_cholesky.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)

# timing tests for increasing matrix sizes for different decomposition algos 
time_chol: $(OBJS)
	$(CPP) $(CPPFLAGS) -o time_chol time_chol.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)

# record matrix reconstruction error for increasing matrix sizes, different algos 
chol_err: $(OBJS)
	$(CPP) $(CPPFLAGS) -o chol_err chol_err.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)

# code for testing applications: Ax=b and Ainverse 
solvers: $(OBJS)
	$(CPP) $(CPPFLAGS) -o solvers solvers.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)

# not needed to run, just testing small matrices again
test_other: $(OBJS)
	$(CPP) $(CPPFLAGS) -o test_other test_other.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)

# rule for generating obj
.cpp.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

all: $(EXECU)

# after tests have been run, move txt files to data directory 
move: 
	mv $(TXT) data/

# remove txt files from data directory 
rmtxt: 
	cd data/ && rm -f $(TXT)

# remove obj files, executables, and txt output files 
clean: 
	rm -f $(OBJS) $(EXECU) $(TXT)