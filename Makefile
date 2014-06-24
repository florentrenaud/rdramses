### Curie ######################################################################
SRCPATH = $(HOME)/bin/src/rdramses.f90
EXEPATH = $(HOME)/bin/rdramses
MPIF = mpif90
F90 = ifort 
FFLAGS = -O3 -cpp
FITS = -L/ccc/cont005/home/gen2192/renaudf/cfitsio/lib/ -lcfitsio
FFT = -L/usr/local/fftw3-3.3.1/lib/ -lfftw3 -lm -I/usr/local/fftw3-3.3.1/include/
################################################################################

### Oxal #######################################################################
#SRCPATH = $(HOME)/bin/src/rdramses.f90
#EXEPATH = $(HOME)/bin/rdramses
#MPIF = mpif90
#F90 = gfortran
#FFLAGS = -O3 -cpp -ffree-line-length-none
#FITS = -L/opt/core-3.1-amd64/cfitsio/3.25/lib -lcfitsio -lm
#FFT = -L/opt/core-3.1-amd64/fftw3/3.2.2/lib/
################################################################################


################################################################################
all: rdramses

rdramses: $(SRCPATH) # full code
	$(F90) $(FFLAGS) $(SRCPATH) -o $(EXEPATH) $(FITS) -DFITS $(FFT) -DFFT

rdramses_nofft: $(SRCPATH) # no FFTW library
	$(F90) $(FFLAGS) $(SRCPATH) -o $(EXEPATH) $(FITS) -DFITS

rdramses_nofits: $(SRCPATH) # no CFITSIO library
	$(F90) $(FFLAGS) $(SRCPATH) -o $(EXEPATH) $(FFT) -DFFT

rdramses_mpi: $(SRCPATH)  # full code, MPI version
	$(MPIF) $(FFLAGS) $(SRCPATH) -o $(EXEPATH)_mpi $(FITS) -DFITS $(FFT) -DFFT -DMPIMODE


rdramses_nothing: $(SRCPATH) # no external library needed
	$(F90) $(FFLAGS) $(SRCPATH) -o $(EXEPATH)


clean :
	rm -f *.o
