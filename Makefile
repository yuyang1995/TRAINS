SHELL=/bin/bash

#CC       = mpiicc -mkl
CC       = mpicc
# OPTIMIZE = -O2 -Wall -finline-functions -fcommon
OPTIMIZE = -O2 -g
# OPTIMIZE += -DDebug

mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))

#------------target system---------
SYSTEM="Linux"

ifeq ($(SYSTEM), "Linux")
NCORE      :=$(grep -c ^processor /proc/cpuinfo)
GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
DNEST_INCL  = -I /home/yuyang/Desktop/Projects/CDNest/
DNEST_LIBS  = -L /home/yuyang/Desktop/Projects/CDNest -ldnest
MPICHINCL     = $(shell pkg-config --cflags mpich) 
MPICHLIB    = $(shell pkg-config --libs mpich)
endif

EXEC     = trains
SRC      = ./src
HEAD	 = ./include
OBJS     = $(SRC)/vars.o $(SRC)/vars_pt.o $(SRC)/main.o $(SRC)/run.o  \
           $(SRC)/read.o $(SRC)/read_pt.o $(SRC)/init.o $(SRC)/range.o \
		   $(SRC)/sim.o $(SRC)/dnest_gen.o  $(SRC)/reconstruct.o \
		   $(SRC)/dnest_pt.o $(SRC)/reconstruct_pt.o $(SRC)/residual.o $(SRC)/MxAvPhase.o\
 		   $(SRC)/dict.o $(SRC)/file.o \
           $(SRC)/system.o $(SRC)/help.o $(SRC)/command_line.o $(SRC)/version.o

INCL =  Makefile $(HEAD)/trains.h $(HEAD)/vars.h $(HEAD)/vars_pt.h \
		$(HEAD)/parset.h $(HEAD)/model.h $(HEAD)/proto.h \
		$(HEAD)/dict.h $(HEAD)/file.h $(HEAD)/system.h $(HEAD)/uthash.h

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(MPICHINCL) $(DNEST_INCL) -I ./include
LIBS     = $(GSL_LIBS) $(MPICHLIB) $(DNEST_LIBS) 

$(EXEC):$(OBJS)
	cd $(SRC)
	$(CC) $(OPTIONS) $(OBJS) $(LIBS) -o $@

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)
