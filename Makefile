goal: makefile.dep
	make bader.x
	make gridchg.x
	make nci.x
	make nat_ibo.x
	make newton.x

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
	-/bin/rm -f *.{o,mod,x,il,a} *__genmod.f90 makefile.dep

makefile.dep: $(shell echo *.f90)
	./make-depend.sh $^ > $@

#
# Explicit dependencies
#
LIBFC += accuracy.o
LIBFC += atoms.o
LIBFC += dgefa.o
LIBFC += dgedi.o
LIBFC += fileio.o
LIBFC += gamess_internal.o
LIBFC += import_gamess.o
LIBFC += lapack.o
LIBFC += lebedev.o
LIBFC += math.o
LIBFC += molden.o
LIBFC += molecular_grid.o
LIBFC += os_integral_operators.o
LIBFC += timer.o
LIBFC += utilities.o

#
# Building the binaries
#
bader.x: bader.o $(LIBFC)
	$(F90) -o bader.x bader.o $(LIBFC) $(LIBS)

gridchg.x: gridchg.o $(LIBFC)
	$(F90) -o gridchg.x gridchg.o $(LIBFC) $(LIBS)

nci.x: nci.o $(LIBFC)
	$(F90) -o nci.x nci.o $(LIBFC) $(LIBS)

nat_ibo.x: nat_ibo.o $(LIBFC)
	$(F90) -o nat_ibo.x nat_ibo.o $(LIBFC) $(LIBS)

newton.x: newton.o $(LIBFC)
	$(F90) -o newton.x newton.o $(LIBFC) $(LIBS)

#
# Automatically-generated dependencies
#
include makefile.dep

