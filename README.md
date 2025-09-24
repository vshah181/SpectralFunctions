# Spectral Function
Plot the surface-pojected spectral functions of wannier90 _hr.dat files using Fortran and python

*Input files:* 
There are four compulsory input files. There is also an optional fifth input file
- First is the master input file.
- Second is the _hr.dat output file from a wannier90 calculation. This needs the filename ’*seedname*_hr.dat'
- Third is the .nnkp file. This needs the filename ’*seedname*.nnkp'
- Fourth is a kpoints file.
- Fifth is a vector potential file, required if constructing a Floquet hamiltonian. The filename is 'vector_potential.dat'
  
*Output files:*
- *seedname*_spec_func.dat. This is the data file with the spectral function for each layer

## Master input file
- This file needs the filename 'INPUT'
1. *seedname* (eg KTaO<sub>3</sub>, SrTiO<sub>3</sub>, BiTeI etc...)
2. The fermi level in electron-volts. This is an optional tag. If specified, this value will be subtracted from the energies when the band structure is plotted 
3. The energy range in which to compute the spectral function
4. The step size in electronvolts for the energy window
5. The broadening factor
6. The bulk switch (1 = do a bulk calculation, 0 = break translational symmetry along a specified direction)
7. The number of layers, to break translational symmetry along an axis. (Ignored if bulk = 1)
8. The direction along which to break symmetry (1 = break transational symmetry along a1 axis, 2 = break translational symmetry along a2 axis. Ignored if bulk = 1)
9. floquet: Whether to do a floquet calculation or not (1 = True, 0 = False)
10. soc: Whether soc is included in the tight-binding hamiltonian (1 = True, 0 = False)
11. basis: [uudd or udud] for the spin. These are different depending on which version of wannier90 was used.
12. figsize (two numbers): the figures size in inches
13. colourmap: which colourmap to plot the spectral function with
14. band_yrange: y limits for the band structure plot
15. bands_plot: whether we need to plot the band structure or not
16. spectra_plot: whether to plot the spectral_function or not
### Example
    seedname           mote2
    e_fermi           -0.5381
    #################################################
    energy_range      -3.3320632588 -0.3320632588
    energy_step        0.0003
    broadening_factor  0.00005
    #################################################
    bulk               1
    nlayers            50
    direction          2
    floquet            1
    #################################################
    soc                .true.
    basis              uudd
    ################# for the plot ##################
    figsize            4 4
    colourmap          inferno
    band_yrange        -0.6 0.6
    bands_plot         1
    spectra_plot       0
*Spaces must be used for separation! Tabs will cause errors!*

## kpoints file
- This file needs the filename 'kpoints'
1. number of k points per path
2. number of high symmetry points
3. high symmetry point 1 (in fractional coordinates) and (optionally) its symbol
4. high symmetry point 2 (in fractional coordinates) and (optionally) its symbol
5. etc...
### Example
    100
    3
    0.10 0.00 0.00 X
    0.00 0.00 0.00 G
    0.10 0.10 0.00 M
In this example, there are 200 kpoints in total and we go along the X - $\Gamma$ - M direction

## vector potential file
- This file needs the filename 'vector_potential.dat'
- For cirulary polarised light we assume a vector potential of the form $(a_0\cos(\omega t), a_0 sin(\omega t+\phi), 0.0)$
- The dimensionless parameter s is used where $s \equiv {e a_0 a \over 2 \hbar c}$
1. s
2. hbar*omega (eV)
3. phi (/pi)
4. maximum order for the floquet expansion
### Example
    s               1.0
    hbar*omega      2.25
    phase_shift/pi  0.0
    max_order       2

# Dependency List
- Fortran compiler (ifort recommended. Also tested on gfortran and cray compiler)
- Blas and lapack implementation (intel math kernel library recommended, also works with cray libsci)
- MPI implementation (mpich or OpenMPI)
- Python $\geq$ 3.7.2
- Numpy
- Matplotlib
