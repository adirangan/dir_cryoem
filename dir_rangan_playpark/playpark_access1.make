
# Compiler
CC=gcc
CCFLAGS_dev= -w -O2 -lpthread -fopenmp -lgslcblas -lfftw3f -lfftw3 -lfinufft -lm -I./dir_h -I./dir_c -L/home/rangan/finufft/lib -I/home/rangan/finufft/include -Wl,-R/home/rangan/finufft/lib -DSINGLE
CCFLAGS=$(CCFLAGS_dev) -D_MONOLITH 

vpath %.c ./dir_c

project=playpark
project_dev=playpark_dev

CLINK=$(CC) -o $(project_dev)

header_dev = ./dir_h/playpark_header.h

sources_dev = array_printf.c \
	array_stats.c \
	boxmuller.c \
	MDA_io.c \
	fillzero.c \
	GLOBAL_pthread.c \
	PNMfile.c \
	updateglobals.c \
	wkspace.c \
	array_transpose.c \
	fftw3_call.c \
	finufft_call.c \
	playpark.c

objects_dev = $(patsubst %.c,./dir_o/%.o,$(sources_dev)) 

all: $(objects_dev) $(project_dev)

$(objects_dev): | ./dir_o

./dir_o:
	@mkdir -p $@

$(project_dev): $(objects_dev) $(header_dev)
	rm -f $(project_dev)
	$(CLINK) $(objects_dev) $(CCFLAGS_dev)

$(project): $(objects_dev) $(header_dev)
	rm -f $(project)
	gcc -w -O2 -lpthread -fopenmp ./dir_c/playpark.c -o $(project) -L./dir_h -L./dir_c -L/home/rangan/finufft/lib -I/home/rangan/finufft/include -I./dir_h -I./dir_c -D_MONOLITH -DSINGLE -lgslcblas -lfftw3f -lfftw3 -lfinufft -lm -Wl,-R/home/rangan/finufft/lib 

./dir_o/%.o : %.c ./dir_h/playpark_header.h
	@echo $< 
	@$(CC) $(CCFLAGS_dev) -c $< -o $@

.c.o: ./dir_h/playpark_header.h
	$(CC) -c $(CCFLAGS_dev) $<

clean: 
	rm -f $(objects_dev)

list_sources: $(sources_dev)
	echo $^

list_objects: $(objects_dev)
	echo $^

