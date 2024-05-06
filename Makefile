CPP = g++
C = gcc
CPPFLAGS = -std=c++20 -I/opt/homebrew/opt/armadillo/include -O1
INCLUDES = -Iinc -I/opt/homebrew/opt/armadillo/include
LDFLAGS = -L/opt/homebrew/opt/armadillo/lib -larmadillo -lblas -llapack
OBJS = cholesky.o matgen.o testing.o
TXT = pivoted_cholesky_times.txt cholesky_times.txt armachol_times.txt cholesky_error.txt pivoted_cholesky_error.txt \
	arma_error.txt LU_times.txt LU_pivot_times.txt cholesky_inv_err.txt piv_cholesky_inv_err.txt \
	LU_decomp_error.txt LU_pivot_decomp_error.txt pivoted_cholesky_times_sing.txt pivoted_cholesky_err_sing.txt \
	cholesky_Axb_err.txt cholesky_piv_Axb_err.txt full_piv_cholesky_error.txt full_pivoted_cholesky_times.txt \
	cholesky_full_piv_Axb_err.txt full_piv_cholesky_inv_err.txt full_pivoted_cholesky_err_sing.txt

EXECU = test_cholesky time_chol chol_err solvers test_other

cholesky.o: cholesky.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c cholesky.cpp -o cholesky.o

matgen.o: matgen.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c matgen.cpp -o matgen.o

testing.o: testing.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c testing.cpp -o testing.o

test_cholesky: $(OBJS)
	$(CPP) $(CPPFLAGS) -o test_cholesky test_cholesky.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)
time_chol: $(OBJS)
	$(CPP) $(CPPFLAGS) -o time_chol time_chol.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)
chol_err: $(OBJS)
	$(CPP) $(CPPFLAGS) -o chol_err chol_err.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)

solvers: $(OBJS)
	$(CPP) $(CPPFLAGS) -o solvers solvers.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)
test_other: $(OBJS)
	$(CPP) $(CPPFLAGS) -o test_other test_other.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)
.cpp.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

all: $(EXECU)
move: 
	mv $(TXT) data/
rmtxt: 
	cd data/ && rm -f $(TXT)
clean: 
	rm -f $(OBJS) $(EXECU) $(TXT)