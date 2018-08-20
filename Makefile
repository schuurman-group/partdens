goal: makefile.dep
	make bader.x
	make ddcharge.x
	make sphatom.x

MAKEFLAGS = -r

.SUFFIXES: .f90 .o .x .c .dep

#
# System-specific overrides
#
include config.mak

#
# Finish the set-up
#
LIBS = $(LAPACK) $(FFTW3)

#
# Compiling and archiving rules
#
.f90.o:
	$(F90) -c $<

dgefa.o:	dgefa.f
	$(F90) -c dgefa.f

dgedi.o:	dgedi.f
	$(F90) -c dgedi.f

.c.o:
	$(CC) -c $<

clean:
	-/bin/rm -f *.{o,mod,x,il,a} checkpoint_{field,main}.* makefile.dep

makefile.dep: $(shell echo *.f90)
	./make-depend.sh $^ > $@

#
# Explicit dependencies
#
LIBFC += lapack.o
LIBFC += math.o
LIBFC += timer.o
LIBFC += dgefa.o
LIBFC += dgedi.o
LIBFC += gamess_internal.o
LIBFC += import_gamess.o
LIBFC += os_integral_operators.o
LIBFC += accuracy.o
LIBFC += atoms.o
LIBFC += lebedev.o
LIBFC += molecular_grid.o
LIBFC += atomdens.o

#
# Building the binaries
#
bader.x: bader.o $(LIBFC)
	$(F90) -o bader.x bader.o $(LIBFC) $(LIBS)

ddcharge.x: ddcharge.o $(LIBFC)
	$(F90) -o ddcharge.x ddcharge.o $(LIBFC) $(LIBS)

sphatom.x: sphatom.o $(LIBFC)
	$(F90) -o sphatom.x sphatom.o $(LIBFC) $(LIBS)

#
# Automatically-generated dependencies
#
include makefile.dep

