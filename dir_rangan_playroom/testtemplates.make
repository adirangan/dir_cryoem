# Barnett adjusted for global include flag options. 6/24/16

# get the global options (barnett 6/23/16)
include make.inc.dist
-include make.inc

PROJECT=int2         # for historical reasons we always like to 
                     # call the executable int2, but you could set it to 
                     # something more descriptive
OBJSUF=o
MODSUF=mod
FLINK=$(FC) -o $(PROJECT)

.f.$(OBJSUF):
	$(FC) -c $(FFLAGS) $<

.f.$(MODSUF):
	$(FC) -c $(FFLAGS) $<

.SUFFIXES: $(MODSUF) .$(OBJSUF) .f .c

vpath %.f .:../spharmproj

FMODS = 

FSRCS = testtemplates.f  \
        template_gen.f \
        template_trans.f \
        gridrouts3d.f \
        spharmroutsall.f \
        spharmroutsnewp.f \
        testfourier.f \
        kspace_model_norm.f \
        testprojeval.f \
        testgreatcircles.f \
        testlsq.f \
        utils.f \
        gaussq.f \
        chebexps.f \
        yrecursion.f \
        lsqsolvespharm.f \
        nufft3df90.f \
        next235.f \
        dfftpack.f \
        ccongrlsq.f \
        prini.f 

#
# object files list
MODS    =  $(FMODS:.f=.$(MODSUF)) 
OBJS    =  $(FMODS:.f=.$(OBJSUF)) $(FSRCS:.f=.$(OBJSUF)) 
#
$(PROJECT):   $(MODS)   $(OBJS)
	rm -f $(PROJECT)
#	$(FLINK) $(OBJS) $(LFLAGS) finufft.a -fopenmp -lstdc++
	$(FLINK) $(OBJS) $(LFLAGS) -L/home/rangan/finufft-master/lib/ -lfinufft -lfftw3 -lfftw3_threads -fopenmp -lstdc++
	./$(PROJECT)
#
clean: 
	rm -f $(OBJS)
# 
list: $(FSRCS)
	echo $^
#
pack: $(FSRCS)
	cat $^ > _tmp_.pack.f
#
