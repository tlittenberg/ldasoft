CC = mpicc

LIBS  = gsl gslcblas m lisa gbmcmc tools omp gsl hdf5

LDFLAGS = -L/Users/tyson/ldasoft/global_fit/lib
INCDIR =/Users/tyson/ldasoft/global_fit/include/

#CCFLAGS += -g -ffast-math -Wall -O3 -ftree-vectorize -std=gnu99 -Werror
CCFLAGS += -g -ffast-math -Wall -O3 -ftree-vectorize -std=gnu99

CCFLAGS += -Xpreprocessor -fopenmp

# include gb_mcmc routines
CCFLAGS += -I$(INCDIR)

OBJS = GalacticBinaryWrapper.o NoiseWrapper.o

all: $(OBJS) global_fit

GalacticBinaryWrapper.o : GalacticBinaryWrapper.c GalacticBinaryWrapper.h
	$(CC) $(CCFLAGS) -c GalacticBinaryWrapper.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

NoiseWrapper.o : NoiseWrapper.c NoiseWrapper.h
	$(CC) $(CCFLAGS) -c NoiseWrapper.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

global_fit : GlobalFit.c $(OBJS)
	$(CC) $(CCFLAGS) -o global_fit GlobalFit.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LDFLAGS) $(LIBS:%=-l%)


clean:
	rm -rf *.dSYM *.o global_fit

