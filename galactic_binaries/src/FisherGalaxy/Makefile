# OpenMP
#CC = gcc-mp-5
CC = gcc
#CCFLAGS = -fopenmp

LIBS  = gsl gslcblas m
#CCFLAGS += -g -ffast-math -Wall -O3 -ftree-vectorize -std=gnu99 -Werror 
CCFLAGS += -g -ffast-math -Wall -O3 -ftree-vectorize -std=gnu99 

# Compile src with git hash
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
CCFLAGS += -DVERSION=\"$(GIT_VERSION)\"

OBJS = arrays.o Subroutines.o

all: $(OBJS) Confusion_Fit Bright_Remove Galaxy OrbitFile Setup Fisher_Galaxy

arrays.o : arrays.c
	$(CC) $(CCFLAGS) -c arrays.c

Subroutines.o : Subroutines.c
	$(CC) $(CCFLAGS) -c Subroutines.c 

Confusion_Fit : Confusion_Fit.c
	$(CC) $(CCFLAGS) -o Confusion_Fit Confusion_Fit.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

Bright_Remove : Bright_Remove.c
	$(CC) $(CCFLAGS) -o Bright_Remove Bright_Remove.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

Galaxy : Galaxy.c
	$(CC) $(CCFLAGS) -o Galaxy Galaxy.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

OrbitFile : OrbitFile.c
	$(CC) $(CCFLAGS) -o OrbitFile OrbitFile.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

Setup : Setup.c
	$(CC) $(CCFLAGS) -o Setup Setup.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)

Fisher_Galaxy : Fisher_Galaxy.c
	$(CC) $(CCFLAGS) -o Fisher_Galaxy Fisher_Galaxy.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)


install : Confusion_Fit Bright_Remove Galaxy OrbitFile Setup Fisher_Galaxy
	install Confusion_Fit Bright_Remove Galaxy OrbitFile Setup Fisher_Galaxy ${HOME}/ldasoft/master/bin/

clean:
	rm -rf *.dSYM *.o Confusion_Fit Bright_Remove Galaxy OrbitFile Setup Fisher_Galaxy
