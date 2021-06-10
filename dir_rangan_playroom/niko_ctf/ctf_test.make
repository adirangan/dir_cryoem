# Marina's general makefile for all tests.
# Alex added global openmp options 6/22/16

# get the global options (barnett 6/23/16)
include make.inc.dist
-include make.inc

PROJECTS=ctf_test_2 

OBJSUF=o
FFLAGS += -g

.f.$(OBJSUF):
	$(FC) -c $(FFLAGS) $<

.SUFFIXES: .$(OBJSUF) .f .c

COMMON_OBJS = \
	ctf_functions.o \

all: $(PROJECTS)


ctf_test: ctf_test.o $(COMMON_OBJS)
	$(FC) -o $@ ctf_test.o $(COMMON_OBJS) $(LFLAGS)

ctf_test_2: ctf_test_2.o $(COMMON_OBJS)
	$(FC) -o $@ ctf_test_2.o $(COMMON_OBJS) $(LFLAGS)

clean: 
	-rm -f *.o ctf_test ctf_test_2

