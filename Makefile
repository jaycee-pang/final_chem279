CPP = g++
C = gcc
CPPFLAGS = -std=c++20 -I/opt/homebrew/opt/armadillo/include
INCLUDES = -Iinc -I/opt/homebrew/opt/eigen/include/eigen3 -I/opt/homebrew/opt/armadillo/include
LDFLAGS = -L/opt/homebrew/opt/armadillo/lib -larmadillo -lblas -llapack
OBJS = cholesky.o matgen.o
TXT = pivoted_cholesky_times.txt cholesky_times.txt armachol_times.txt cholesky_error.txt pivoted_cholesky_error.txt \
	arma_error.txt

EXECU = test_cholesky time_chol chol_err solvers

cholesky.o: cholesky.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c cholesky.cpp -o cholesky.o

matgen.o: matgen.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c matgen.cpp -o matgen.o

test_cholesky: $(OBJS)
	$(CPP) $(CPPFLAGS) -o test_cholesky test_cholesky.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)
time_chol: $(OBJS)
	$(CPP) $(CPPFLAGS) -o time_chol time_chol.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)
chol_err: $(OBJS)
	$(CPP) $(CPPFLAGS) -o chol_err chol_err.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)

solvers: $(OBJS)
	$(CPP) $(CPPFLAGS) -o solvers solvers.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)

.cpp.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

all: $(EXECU)
move: 
	mv $(TXT) data/
rmtxt: 
	cd data/ && rm -f $(TXT)
clean: 
	rm -f $(OBJS) $(EXECU) $(TXT)