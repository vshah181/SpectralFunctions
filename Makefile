PROGRAM = spectral_function
FTN = mpif90
FLIB = -lblas -llapack
FFLAGS = -O4 -Wall -Wextra

CONSTNT = constants.mod
FILE_IO = file_parsing.mod
PEIERLS = peierls_substitution.mod
MODULES = constants.mod file_parsing.mod peierls_substitution.mod
OBJMDLS = constants.o file_parsing.o peierls_substitution.o

OSUBRTN = add_potential.o energy_window.o floquet_expansion.o fourier_coefficient.o fourier_transform.o greens_function.o kpath_maker.o main.o
SSUBRTN = add_potential.f90 energy_window.f90 floquet_expansion.f90 fourier_coefficient.f90 fourier_transform.f90 greens_function.f90 kpath_maker.f90 main.f90

OBJFILE = $(OBJMDLS) $(OSUBRTN)

all: $(PROGRAM)

# Build the program
$(PROGRAM): $(OBJFILE) $(OSUBRTN) $(MODULES)
	$(FTN) $(OBJFILE) -o $@ $(FLIB) $(FFLAGS)

# Make the subroutines (they need the modules)
$(OSUBRTN): %.o: %.f90 $(MODULES)
	$(FTN) -c $< -o $@ $(FLIB) $(FFLAGS)

# Make all the modules
$(OBJMDLS): $(PEIERLS) $(FILE_IO) $(CONSTNT);
$(PEIERLS): peierls_substitution.f90 $(CONSTNT) $(FILE_IO)
	$(FTN) -c $< $(FLIB) $(FFLAGS)
$(FILE_IO): file_parsing.f90 $(CONSTNT)
	$(FTN) -c $< $(FLIB) $(FFLAGS)
$(CONSTNT): constants.f90
	$(FTN) -c $< $(FLIB) $(FFLAGS)


clean:
	rm $(OBJFILE) $(PROGRAM) $(MODULES)
