# OpenMP
#CC = gcc-mp-5
CC = clang-mp-9.0
CCFLAGS = -fopenmp 

LIBS  = gsl gslcblas m fftw3
CCFLAGS += -O3 -g -ffast-math -Wall -ftree-vectorize -std=gnu99 -Werror

# Compile src with git hash
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
CCFLAGS += -DVERSION=\"$(GIT_VERSION)\"

OBJS = LISA.o GalacticBinaryIO.o GalacticBinaryModel.o GalacticBinaryWaveform.o GalacticBinaryMath.o GalacticBinaryData.o GalacticBinaryPrior.o GalacticBinaryProposal.o GalacticBinaryFStatistic.o 

all: gitversion.h $(OBJS) gb_mcmc gb_match gb_catalog gb_mcmc_chirpmass gb_mcmc_brans_dicke gb_residual libgbmcmc.a

gitversion.h: ../../.git/HEAD ../../.git/index
	echo "const char *gitversion = \"$(shell git rev-parse HEAD)\";" > $@

LISA.o : LISA.c LISA.h
	$(CC) $(CCFLAGS) -c LISA.c 

GalacticBinaryIO.o : GalacticBinaryIO.c GalacticBinaryIO.h LISA.h
	$(CC) $(CCFLAGS) -c GalacticBinaryIO.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

GalacticBinaryMath.o : GalacticBinaryMath.c GalacticBinaryMath.h GalacticBinary.h
	$(CC) $(CCFLAGS) -c GalacticBinaryMath.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

GalacticBinaryPrior.o : GalacticBinaryPrior.c GalacticBinaryPrior.h GalacticBinary.h
	$(CC) $(CCFLAGS) -c GalacticBinaryPrior.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

GalacticBinaryData.o : GalacticBinaryData.c GalacticBinaryData.h GalacticBinary.h GalacticBinaryModel.o GalacticBinaryMath.o GalacticBinaryIO.o LISA.o 
	 $(CC) $(CCFLAGS) -c GalacticBinaryData.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

GalacticBinaryWaveform.o : GalacticBinaryWaveform.c GalacticBinaryWaveform.h GalacticBinaryMath.o Constants.h LISA.o 
	$(CC) $(CCFLAGS) -c GalacticBinaryWaveform.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

GalacticBinaryFStatistic.o : GalacticBinaryFStatistic.c GalacticBinaryWaveform.c GalacticBinary.h Constants.h LISA.o
	$(CC) $(CCFLAGS) -c GalacticBinaryFStatistic.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

GalacticBinaryModel.o : GalacticBinaryModel.c GalacticBinaryModel.h LISA.o GalacticBinary.h GalacticBinaryMath.o GalacticBinaryWaveform.o
	$(CC) $(CCFLAGS) -c GalacticBinaryModel.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

GalacticBinaryProposal.o : GalacticBinaryProposal.c GalacticBinary.h GalacticBinaryWaveform.o GalacticBinaryIO.o GalacticBinaryModel.o Constants.h
	$(CC) $(CCFLAGS) -c GalacticBinaryProposal.c $(INCDIR:%=-I%) $(LIBDIR:%=-L%)

gb_mcmc : GalacticBinaryMCMC.c $(OBJS) GalacticBinary.h
	$(CC) $(CCFLAGS) -o gb_mcmc GalacticBinaryMCMC.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

gb_catalog : GalacticBinaryCatalog.c $(OBJS) GalacticBinary.h
	$(CC) $(CCFLAGS) -o gb_catalog GalacticBinaryCatalog.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

gb_residual : GalacticBinaryResidual.c $(OBJS) GalacticBinary.h
	$(CC) $(CCFLAGS) -o gb_residual GalacticBinaryResidual.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

gb_mcmc_chirpmass : gb_mcmc_chirpmass.c 
	$(CC) -o gb_mcmc_chirpmass gb_mcmc_chirpmass.c $(LIBS:%=-l%)

gb_match : GalacticBinaryMatch.c $(OBJS) GalacticBinary.h
	$(CC) $(CCFLAGS) -o gb_match GalacticBinaryMatch.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

gb_mcmc_brans_dicke : gb_mcmc_brans_dicke.c
	$(CC) -o gb_mcmc_brans_dicke gb_mcmc_brans_dicke.c $(LIBS:%=-l%)

libgbmcmc.a : $(OBJS)
#	$(CC) -shared -o libgb.so $(OBJS) $(LIBS:%=-l%)
	ar rcs libgbmcmc.a $(OBJS)

install : gitversion.h gb_mcmc gb_match gb_catalog gb_mcmc_chirpmass gb_mcmc_brans_dicke gb_residual
	install gb_mcmc gb_match gb_catalog gb_mcmc_chirpmass gb_mcmc_brans_dicke ${HOME}/ldasoft/master/bin/

clean:
	rm *.o *.a gb_mcmc gb_match gb_catalog gb.so gb_mcmc_chirpmass gb_mcmc_brans_dicke gb_residual gitversion.h

