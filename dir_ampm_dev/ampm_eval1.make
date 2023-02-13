
# Compiler
CC=gcc
CCFLAGS_ori= -fPIC -w -O2 -D_WITHOUT_MEX -fopenmp -lm -lgsl -I./dir_h -I/usr/local/include
CCFLAGS_intrin= -D_AVX -mavx -D_FMA -mfma
CCFLAGS_fftw3= -lfftw3_omp -lfftw3 -lfftw3f_omp -lfftw3f 
CCFLAGS_cblas= -D_CBLAS -Wl,-R -Wl,/usr/lib/x86_64-linux-gnu -lgslcblas
CCFLAGS_complex= -D_COMPLEX
CCFLAGS_dev = $(CCFLAGS_ori) $(CCFLAGS_intrin) $(CCFLAGS_fftw3) $(CCFLAGS_cblas) $(CCFLAGS_complex)
CCFLAGS=$(CCFLAGS_dev) -D_MONOLITH

vpath %.c ./dir_c

project=ampm
project_dev=ampm_dev
project_lib=ampm_lib
project_ar=libampm

CLINK_dev=$(CC) -o $(project_dev)
CLINK=$(CC) -o $(project)

header_dev = ./dir_h/ampm_header.h

sources_dev = ping_call.c \
	malloc1_call.c \
	update_global_call.c \
	array_printf_call.c \
	array_stats_call.c \
	array_extract_call.c \
	array_transpose_call.c \
	rand_call.c \
	MDA_io_call.c \
	local_tictoc_call.c \
	GLOBAL_tictoc_call.c \
	complex_interleave_call.c \
	transpose_call.c \
	matrix_multiply_call.c \
	gsl_legpts_call.c \
	gsl_sf_legendre_sphPlm_call.c \
	sample_sphere_call.c \
	mex_ampmh_X_wSM_call.c \
	main_driver.c \

objects_dev = $(patsubst %.c,./dir_o/%.o,$(sources_dev)) 

all: $(objects_dev) $(project_dev)

$(objects_dev): | ./dir_o

./dir_o:
	@mkdir -p $@

$(project_dev): $(objects_dev) $(header_dev)
	rm -f $(project_dev)
	$(CLINK_dev) $(objects_dev) $(CCFLAGS_dev)

$(project): $(objects_dev) $(header_dev)
	rm -f $(project)
	$(CLINK) $(objects_dev) $(CCFLAGS)

$(project_lib): $(objects_dev) $(header_dev)
	rm -f $(project_lib).so
	gcc -shared -o $(project_lib).so $(objects_dev) $(CCFLAGS_dev)

$(project_ar): $(objects_dev) $(header_dev)
	rm -f $(project_ar).a
	ar -rcs $(project_ar).a $(objects_dev)

./dir_o/%.o : %.c ./dir_h/ampm_header.h
	@echo $< 
	@$(CC) $(CCFLAGS_dev) -c $< -o $@

.c.o: ./dir_h/ampm_header.h
	$(CC) -c $(CCFLAGS_dev) $<

clean: 
	rm -f $(objects_dev)

list_sources: $(sources_dev)
	echo $^

list_objects: $(objects_dev)
	echo $^

