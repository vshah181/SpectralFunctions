OBJFILES = constants.o file_parsing.o fourier_transform.o kpath_maker.o greens_function.o main.o energy_window.o add_potential.o
PROGRAM = spectral_function
FTN = mpif90
FFLAGS = -O4 -llapack -lblas -ld_classic

all: $(PROGRAM)

$(PROGRAM): $(OBJFILES)
	$(FTN) -o $(PROGRAM) $(OBJFILES) $(FFLAGS)

%.o: %.f90
	$(FTN) -c $< -o $@ $(FFLAGS)

clean:
	rm $(OBJFILES) $(PROGRAM)
