#-----------------------------------------------------------------------
# Compiler flags
#-----------------------------------------------------------------------

#
# gfortran
#
F90	= gfortran
F77	= gfortran
CC	= gcc
F90OPTS = -cpp -g -ffixed-line-length-none -ffree-line-length-none -fopenmp -O3 -fbacktrace
CCOPTS  = -g -O0

# External libraries
LIBS= -lblas -llapack

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

BLOCKDIAG = blockdiag/bdglobal.o \
	blockdiag/adt.o \
	blockdiag/overlaps.o \
	blockdiag/blockdiag.o

DISPLACE = displace/dispglobal.o \
	displace/displace.o

KDC = kdc/kdcglobal.o \
	kdc/opermod.o \
	kdc/parinfo.o \
	kdc/kdc.o

OBJECTS_BLOCKDIAG = $(MULTI) \
	$(INCLUDE) \
	$(IOMODULES) \
	$(UTILITIES) \
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
	timingmod.o \
	iomod.o \
	parsemod.o \
	utils.o \
	bdglobal.o \
	adt.o \
	overlaps.o \
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
	$(KDC)
OBJ_KDC = constants.o \
	channels.o \
	sysinfo.o \
	iomod.o \
	parsemod.o \
	timingmod.o \
	utils.o \
	ioqc.o \
	symmetry.o \
	kdcglobal.o \
	opermod.o \
	parinfo.o \
	kdc.o

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

%.o: %.f90
	$(F90) -c $(F90OPTS) $<

%.o: %.f
	$(F77) -c $(F77OPTS) $<

%.o: %.c
	$(CC) $(CCOPTS)  -c $<

clean_all:
	rm -f *.o *~ *.mod
