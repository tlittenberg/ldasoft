# OpenMP
#CC = gcc-mp-5
CC = gcc
#CCFLAGS = -fopenmp

LIBS  = gsl gslcblas m gbmcmc # fftw3
LDFLAGS = -L../../gbmcmc/src
#CCFLAGS += -g -ffast-math -Wall -O3 -ftree-vectorize -std=gnu99 -Werror
CCFLAGS += -g -ffast-math -Wall -O3 -ftree-vectorize -std=gnu99

# Compile src with git hash
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
CCFLAGS += -DVERSION=\"$(GIT_VERSION)\"

# include gb_mcmc routines
CCFLAGS += -I../../gbmcmc/src/

OBJS = arrays.o Subroutines.o

all: $(OBJS) Confusion_Fit Bright_Remove Galaxy OrbitFile Setup Fisher_Galaxy Full_Residual

arrays.o : arrays.c
	$(CC) $(CCFLAGS) -c arrays.c

Subroutines.o : Subroutines.c Detector.h
	$(CC) $(CCFLAGS) -c Subroutines.c

Confusion_Fit : Confusion_Fit.c Subroutines.o arrays.o
	$(CC) $(CCFLAGS) -o Confusion_Fit Confusion_Fit.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LDFLAGS) $(LIBS:%=-l%)

Full_Residual : Full_Residual.c Subroutines.o arrays.o
	$(CC) $(CCFLAGS) -o Full_Residual Full_Residual.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LDFLAGS) $(LIBS:%=-l%)

Bright_Remove : Bright_Remove.c Subroutines.o arrays.o
	$(CC) $(CCFLAGS) -o Bright_Remove Bright_Remove.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LDFLAGS) $(LIBS:%=-l%)

Galaxy : Galaxy.c Subroutines.o arrays.o
	$(CC) $(CCFLAGS) -o Galaxy Galaxy.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LDFLAGS) $(LIBS:%=-l%)

OrbitFile : OrbitFile.c Subroutines.o arrays.o
	$(CC) $(CCFLAGS) -o OrbitFile OrbitFile.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LDFLAGS) $(LIBS:%=-l%)

Setup : Setup.c Detector.h
	$(CC) $(CCFLAGS) -o Setup Setup.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LDFLAGS) $(LIBS:%=-l%)

Fisher_Galaxy : Fisher_Galaxy.c Subroutines.o arrays.o
	$(CC) $(CCFLAGS) -o Fisher_Galaxy Fisher_Galaxy.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LDFLAGS) $(LIBS:%=-l%)


install : Confusion_Fit Bright_Remove Galaxy OrbitFile Setup Fisher_Galaxy Full_Residual
	install Confusion_Fit Bright_Remove Galaxy OrbitFile Setup Fisher_Galaxy Full_Residual ${HOME}/ldasoft/master/bin/

clean:
	rm -rf *.dSYM *.o Confusion_Fit Bright_Remove Galaxy OrbitFile Setup Fisher_Galaxy Full_Residual
