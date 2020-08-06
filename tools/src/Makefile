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


all: $(OBJS) gaussian_mixture_model

gaussian_mixture_model : gaussian_mixture_model.c
	$(CC) $(CCFLAGS) -o gaussian_mixture_model gaussian_mixture_model.c $(LIBS:%=-l%)


install : gaussian_mixture_model
	install gaussian_mixture_model ${HOME}/ldasoft/master/bin/

clean:
	rm -rf *.dSYM *.o gaussian_mixture_model
