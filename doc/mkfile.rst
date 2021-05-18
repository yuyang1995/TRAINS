****************
Makefile
****************

Makefile is used to compile and link the code. 
There are some important configurations worth mentioning::

  SYSTEM="Linux"

  ifeq ($(SYSTEM), "Linux")
  GSL_INCL    = $(shell pkg-config --cflags gsl) 
  GSL_LIBS    = $(shell pkg-config --libs gsl)

  DNEST_INCL  = -I /home/yuyang/Desktop/Projects/CDNest/
  DNEST_LIBS  = -L /home/yuyang/Desktop/Projects/CDNest -ldnest

  MPICHINCL   = $(shell pkg-config --cflags mpich) 
  MPICHLIB    = $(shell pkg-config --libs mpich)
  endif

Locations of header and library files of libraries should be adjusted according to where they are installed.