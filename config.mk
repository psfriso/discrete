DEBUG =
#DEBUG = -g -fbacktrace -fbounds-check 
#DEBUG = -g -fbacktrace -fbounds-check -DDEBUG
DEBUGNOCHECK = 
#DEBUGNOCHECK = -g -fbacktrace -DDEBUG
PROFILE =
#PROFILE = -pg
#PARALLEL =
PARALLEL = -fopenmp -DOMP
#
OPTIONS = -O3 -fimplicit-none -funroll-loops -mtune=generic -Wall -ftree-vectorizer-verbose=1 -ffast-math  
LOPTIONS= -static-libgfortran  
F77= gfortran
CPP= cpp
LIBDIR = ../lib
CFLAGS = -c $(DEBUG) $(PROFILE) -I$(LIBDIR) $(PARALLEL)
CFLAGSNODEBUG = -c $(PROFILE) -I$(LIBDIR) $(PARALLEL)
CFLAGSNOCHECK = -c $(DEBUGNOCHECK) $(PROFILE) -I$(LIBDIR) $(PARALLEL)
LFLAGS = $(PROFILE) -L$(LIBDIR) $(PARALLEL)
ENDFLAG = -fconvert=little-endian
