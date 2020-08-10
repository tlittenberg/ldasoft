CC = gcc

LIBS  = gsl gslcblas m 
#LDFLAGS = -L../../gbmcmc/src
#CCFLAGS += -g -ffast-math -Wall -O3 -ftree-vectorize -std=gnu99 -Werror
CCFLAGS += -g -ffast-math -Wall -O3 -ftree-vectorize -std=gnu99

# Compile src with git hash
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
CCFLAGS += -DVERSION=\"$(GIT_VERSION)\"

# include gb_mcmc routines
#CCFLAGS += -I../../gbmcmc/src/

OBJS = GMM_with_EM.o

all: $(OBJS) gaussian_mixture_model libtools.a

GMM_with_EM.o : GMM_with_EM.c GMM_with_EM.h
	$(CC) $(CCFLAGS) -c GMM_with_EM.c

gaussian_mixture_model : gaussian_mixture_model.c
	$(CC) $(CCFLAGS) -o gaussian_mixture_model gaussian_mixture_model.c $(OBJS) $(LIBS:%=-l%)

libtools.a : $(OBJS)
	ar rcs libtools.a $(OBJS)

install : gaussian_mixture_model
	install gaussian_mixture_model ${HOME}/ldasoft/master/bin/

clean:
	rm -rf *.dSYM *.o gaussian_mixture_model
