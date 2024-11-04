#-----------------------------------------------------------------------
# Compiler flags
#-----------------------------------------------------------------------

#
# gfortran
#
#F90	= gfortran
#F77	= gfortran
#CC	= gcc
#F90OPTS = -cpp -g -ffixed-line-length-none -ffree-line-length-none -fopenmp -O3 -fbacktrace
#CCOPTS  = -g -O0

#
# ifort
#
F90	 = ifort
F77	 = ifort
CC	 = icc
F90OPTS = -cpp -g -free -fopenmp -traceback -O3 -diag-disable 8290 -diag-disable 8291 -diag-disable 10448
CCOPTS  = -g -O0

# External libraries
LIBS= -lblas -llapack
#LIBS = -mkl

#-----------------------------------------------------------------------
# Define object files
#-----------------------------------------------------------------------
MULTI = multi/accuracy.o \
	multi/printing.o \
	multi/timer.o \
	multi/lapack.o \
	multi/dgefa.o \
	multi/dgedi.o \
	multi/math.o \
	multi/os_integral_operators.o \
	multi/gamess_internal.o \
	multi/import_gamess.o 

INCLUDE = include/constants.o \
	include/channels.o \
	include/sysinfo.o

IOMODULES=iomodules/iomod.o \
	iomodules/parsemod.o \

IOQC = ioqc/ioqc.o

UTILITIES=utilities/timingmod.o \
	utilities/utils.o

SYMMETRY=symmetry/symmetry.o

POTFUNCS=potfuncs/parameters.o \
	potfuncs/potfuncs.o

OPT=opt/opt.o

BLOCKDIAG = blockdiag/bdglobal.o \
	blockdiag/adt.o \
	blockdiag/mooverlaps.o \
	blockdiag/wfoverlaps.o \
	blockdiag/tamura.o \
	blockdiag/pacher.o \
	blockdiag/detparsing.o \
	blockdiag/blockdiag.o

DISPLACE = displace/dispglobal.o \
	displace/displace.o

KDC = kdc/kdcglobal.o \
	kdc/cartgrad.o \
	kdc/opermod.o \
	kdc/parinfo.o \
	kdc/nmeqmod.o \
	kdc/transform.o \
	kdc/kdc.o

PLTKDC = pltkdc/pltglobal.o \
	pltkdc/pltkdc.o

MERGEOP = mergeop/mergeglobal.o \
	kdc/kdcglobal.o \
	kdc/opermod.o \
	mergeop/mergeop.o

OBJECTS_BLOCKDIAG = $(MULTI) \
	$(INCLUDE) \
	$(IOMODULES) \
	$(UTILITIES) \
	$(IOQC) \
	$(BLOCKDIAG)
OBJ_BLOCKDIAG = accuracy.o \
	printing.o \
	timer.o \
	lapack.o \
	dgefa.o \
	dgedi.o \
	math.o \
	os_integral_operators.o \
	gamess_internal.o \
	import_gamess.o \
	constants.o \
	channels.o \
	sysinfo.o \
	timingmod.o \
	iomod.o \
	parsemod.o \
	utils.o \
	ioqc.o \
	bdglobal.o \
	adt.o \
	mooverlaps.o \
	wfoverlaps.o \
	tamura.o \
	pacher.o \
	detparsing.o \
	blockdiag.o

OBJECTS_DISPLACE = $(INCLUDE) \
	$(IOMODULES) \
	$(UTILITIES) \
	$(IOQC) \
	$(SYMMETRY) \
	$(DISPLACE)
OBJ_DISPLACE = constants.o \
	channels.o \
	sysinfo.o \
	iomod.o \
	parsemod.o \
	timingmod.o \
	utils.o \
	ioqc.o \
	symmetry.o \
	dispglobal.o \
	displace.o

OBJECTS_KDC = $(INCLUDE) \
	$(IOMODULES) \
	$(UTILITIES) \
	$(IOQC) \
	$(SYMMETRY) \
	$(POTFUNCS) \
	$(OPT) \
	$(KDC)
OBJ_KDC = constants.o \
	channels.o \
	parameters.o \
	sysinfo.o \
	iomod.o \
	parsemod.o \
	timingmod.o \
	utils.o \
	ioqc.o \
	symmetry.o \
	potfuncs.o \
	opt.o \
	kdcglobal.o \
	cartgrad.o \
	opermod.o \
	parinfo.o \
	nmeqmod.o \
	transform.o \
	kdc.o

OBJECTS_PLTKDC = $(INCLUDE) \
	$(IOMODULES) \
	$(SYMMETRY) \
	$(POTFUNCS) \
	$(PLTKDC)

OBJ_PLTKDC = constants.o \
	channels.o \
	parameters.o \
	sysinfo.o \
	iomod.o \
	parsemod.o \
	symmetry.o \
	potfuncs.o \
	pltglobal.o \
	pltkdc.o

OBJECTS_MERGEOP = $(INCLUDE) \
	$(IOMODULES) \
	$(SYMMETRY) \
	$(POTFUNCS) \
	$(MERGEOP)

OBJ_MERGEOP = constants.o \
	channels.o \
	parameters.o \
	symmetry.o \
	sysinfo.o \
	iomod.o \
	parsemod.o \
	mergeglobal.o \
	kdcglobal.o \
	opermod.o \
	mergeop.o

#-----------------------------------------------------------------------
# Rules to create the program
#-----------------------------------------------------------------------
blockdiag: $(OBJECTS_BLOCKDIAG)
	$(F90) $(F90OPTS) $(OBJ_BLOCKDIAG) $(LIBS) -o bin/blockdiag.x
	rm -f *.o *~ *.mod 2>/dev/null

displace: $(OBJECTS_DISPLACE)
	$(F90) $(F90OPTS) $(OBJ_DISPLACE) $(LIBS) -o bin/displace.x
	rm -f *.o *~ *.mod 2>/dev/null

kdc: $(OBJECTS_KDC)
	$(F90) $(F90OPTS) $(OBJ_KDC) $(LIBS) -o bin/kdc.x
	rm -f *.o *~ *.mod 2>/dev/null

pltkdc: $(OBJECTS_PLTKDC)
	$(F90) $(F90OPTS) $(OBJ_PLTKDC) $(LIBS) -o bin/pltkdc.x
	rm -f *.o *~ *.mod 2>/dev/null

mergeop: $(OBJECTS_MERGEOP)
	$(F90) $(F90OPTS) $(OBJ_MERGEOP) $(LIBS) -o bin/mergeop.x
	rm -f *.o *~ *.mod 2>/dev/null

%.o: %.f90
	$(F90) -c $(F90OPTS) $<

%.o: %.f
	$(F77) -c $(F77OPTS) $<

%.o: %.c
	$(CC) $(CCOPTS)  -c $<

clean_all:
	rm -f *.o *~ *.mod
