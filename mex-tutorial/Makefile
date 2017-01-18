# Linux.
MEXSUFFIX  = mexglx
MATLABHOME = /cs/local/generic/lib/pkg/matlab-7.3
MEX        = mex
CXX        = g++-3.4.5
CFLAGS     = -fPIC -ansi -pthread -DMX_COMPAT_32 -DMATLAB_MEX_FILE
LIBPATHS   = -L$(HOME)/lib/gcc-3.4.5

# Mac OS X.
# MEXSUFFIX  = mexmac
# MATLABHOME = /Applications/MATLAB72
# MEX        = mex
# CXX        = g++
# CFLAGS     = -fno-common -no-cpp-precomp -fexceptions 

TARGETS   = normpdf.$(MEXSUFFIX)
LIBS      = -lgsl -lgslcblas -lm
MEXFLAGS  = -cxx CC='$(CXX)' CXX='$(CXX)' LD='$(CXX)'
INCLUDE   = -I$(HOME)/include -I$(MATLABHOME)/extern/include

all: $(TARGETS)

normpdf.o: normpdfgsl.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -c -o normpdf.o $^

normpdf.$(MEXSUFFIX): normpdf.o
	$(MEX) $(MEXFLAGS) $(LIBPATHS) $(LIBS) -output normpdf $^

clean:
	rm -f $(TARGETS) normpdf.o

