# OpenMP
CC = gcc-mp-5
CCFLAGS = -fopenmp

LIBS  = gsl gslcblas m
CCFLAGS += -g -Wall -std=gnu99

# Compile src with git hash
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
CCFLAGS += -DVERSION=\"$(GIT_VERSION)\"

OBJS = LISA.o GalacticBinaryIO.o GalacticBinaryModel.o GalacticBinaryWaveform.o GalacticBinaryMath.o GalacticBinaryData.o

all: $(OBJS) gb_mcmc

LISA.o : LISA.c LISA.h
	$(CC) $(CCFLAGS) -c LISA.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

GalacticBinaryIO.o : GalacticBinaryIO.c GalacticBinaryIO.h LISA.h
	$(CC) $(CCFLAGS) -c GalacticBinaryIO.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

GalacticBinaryMath.o : GalacticBinaryMath.c GalacticBinaryMath.h GalacticBinary.h
	$(CC) $(CCFLAGS) -c GalacticBinaryMath.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

GalacticBinaryData.o : GalacticBinaryData.c GalacticBinaryData.h GalacticBinary.h GalacticBinaryModel.o GalacticBinaryWaveform.o GalacticBinaryMath.o LISA.o
	 $(CC) $(CCFLAGS) -c GalacticBinaryData.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

GalacticBinaryWaveform.o : GalacticBinaryWaveform.c GalacticBinaryWaveform.h GalacticBinaryMath.o Constants.h LISA.o GalacticBinaryModel.o
	$(CC) $(CCFLAGS) -c GalacticBinaryWaveform.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

GalacticBinaryModel.o : GalacticBinaryModel.c GalacticBinaryModel.h
	$(CC) $(CCFLAGS) -c GalacticBinaryModel.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

gb_mcmc : GalacticBinaryMCMC.c $(OBJS)
	$(CC) $(CCFLAGS) -o gb_mcmc GalacticBinaryMCMC.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

clean:
	rm *.o gb_mcmc
