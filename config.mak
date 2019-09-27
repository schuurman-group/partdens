SHELL = bash
CC  = gcc -Wall -m64 -march=k8 -msse3 -mfpmath=sse -O3
F90 = ifort -warn -nofor_main -ipo -O3 -no-prec-div -qopenmp -assume cc_omp -complex_limited_range \
            -debug extended -traceback
F90L = $(F90)
LAPACK = -mkl=sequential
FFTW3 =  -L/usr/local/fftw/3.3.3/intel/impi/lib64 -lfftw3 -lfftw3f

fftw.o:	fftw_fftw3.f90 accuracy.o
	$(F90) -c fftw_fftw3.f90 -o fftw.o
