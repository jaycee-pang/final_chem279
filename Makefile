CPP = g++
C = gcc
CPPFLAGS = -std=c++20 -I/opt/homebrew/opt/armadillo/include
INCLUDES = -Iinc -I/opt/homebrew/opt/eigen/include/eigen3 -I/opt/homebrew/opt/armadillo/include
LDFLAGS = -L/opt/homebrew/opt/armadillo/lib -larmadillo -lblas -llapack
OBJS = cholesky.o matgen.o


EXECU = test_cholesky

cholesky.o: cholesky.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c cholesky.cpp -o cholesky.o

matgen.o: matgen.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c matgen.cpp -o matgen.o

test_cholesky: $(OBJS)
	$(CPP) $(CPPFLAGS) -o test_cholesky test_cholesky.cpp $(OBJS) $(INCLUDES) $(LDFLAGS)

.cpp.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

all: $(EXECU)

clean: 
	rm -f $(OBJS) $(EXECU)