
# Compiler
CC=gcc
CCFLAGS_dev= -fPIC -w -O2 -mavx -mfma -lpthread -fopenmp -lgslcblas -lfftw3f -lfftw3 -lfinufft -lm -I./dir_h -I./dir_c -I/home/rangan/dir_finufft/include -L/home/rangan/dir_finufft/lib  -Wl,-R -Wl,/home/rangan/dir_finufft/lib -Wl,-R -Wl,/usr/lib/x86_64-linux-gnu -DSINGLE
CCFLAGS=$(CCFLAGS_dev) -D_MONOLITH 

vpath %.c ./dir_c

project=playpark
project_dev=playpark_dev
project_lib=playpark_lib

CLINK=$(CC) -o $(project_dev)

header_dev = ./dir_h/playpark_header.h

sources_dev = array_printf_call.c \
	array_stats_call.c \
	boxmuller.c \
	MDA_io_call.c \
	fillzero.c \
	GLOBAL_pthread.c \
	PNMfile.c \
	updateglobals.c \
	wkspace.c \
	array_transpose_call.c \
	array_extract_call.c \
	array_mean_center_call.c \
	fftw3_call.c \
	finufft_call.c \
	dp_ps_call.c \
	get_xdrop_logscale_call.c \
	quicksort_call.c \
	dexcluster_nonbinary_call.c \
	nelder_mead_call.c \
	gumbel_call.c \
	array_orth_call.c \
	erfcln_call.c \
	find_internal_maximum_call.c \
	dexcluster_nonbinary_recursive_call.c \
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
	gcc -fPIC -w -O2 -mavx -mfma -lpthread -fopenmp ./dir_c/playpark.c -o $(project) -L./dir_h -L./dir_c -L/home/rangan/dir_finufft/lib -I/home/rangan/dir_finufft/include -I./dir_h -I./dir_c -D_MONOLITH -DSINGLE -lgslcblas -lfftw3f -lfftw3 -lfinufft -lm -Wl,-R -Wl,/home/rangan/dir_finufft/lib -Wl,-R -Wl,/usr/lib/x86_64-linux-gnu 

$(project_lib): $(objects_dev) $(header_dev)
	rm -f $(project_lib).so
	gcc -shared -fPIC -w -O2 -mavx -mfma -lpthread -fopenmp -o $(project_lib).so $(objects_dev) -L./dir_h -L./dir_c -L/home/rangan/dir_finufft/lib -I/home/rangan/dir_finufft/include -I./dir_h -I./dir_c -D_MONOLITH -DSINGLE -lgslcblas -lfftw3f -lfftw3 -lfinufft -lm -Wl,-R -Wl,/home/rangan/dir_finufft/lib -Wl,-R -Wl,/usr/lib/x86_64-linux-gnu 

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

