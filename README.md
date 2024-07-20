# Spectral Function
Plot the surface-pojected spectral functions of wannier90 _hr.dat files using Fortran and gnuplot

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
2. The number of layers, to break translational symmetry along an axis.
3. The direction along which to break symmetry (1=x, 2=y, 3=z)
4. The basis, the current version of wannier90 the hr file in the up, down, up down... basis whereas the old version writes up, up ..., up, down, down, ..., down
5. The fermi level in electron-volts. This is an optional tag. If specified, this value will be subtracted from the energies when the band structure is plotted 
6. The energy range in which to compute the spectral function
7. The broadening factor
8. The step size in electronvolts for the energy window
9. The output figure size in inches (this will only be read by the python script)
10. Whether or not to do construct a Floquet-Hamilotonian. This requires an addition input file: vector_potential.dat
### Example
    seedname           EuPb
    nlayers            30
    basis              uudd
    e_fermi            6.6393567
    energy_range       5.6  8.0
    broadening_factor  0.003
    energy_step        0.0003
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
- The dimensionless parameter s is used where $s \equiv {e a_0 a \ over 2 \hbar c}$
1. s
2. hbar*omega (eV)
3. phi (/pi)
4. maximum order for the floquet expansion
### Example
    s               1.0
    hbar*omega      2.25
    phase_shift/pi  0.0
    max_order       2

### NOTE: 
The python script requires numpy and matplotlib to work properly.
