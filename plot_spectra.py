import sys
import numpy as np
import matplotlib.pyplot as plt

#TODO This whole script would be better if I made it object-oriented 
plt.rcParams["font.family"] = "Arial"
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Arial"
plt.rcParams["mathtext.it"] = "Arial:italic"
plt.rcParams["mathtext.bf"] = "Arial:bold"
plt.rcParams['svg.fonttype'] = 'none'
def read_layer_index(nlayers, min_layers):
    """
    Read the required layers from the user input.
    :param nlayers: Number of layers in the modelled system (integer).
    :param min_layers: Minimum acceptable layer index (integer).
    :return: layer_index: the layer we need to draw a plot for (integer).
    """
    user_input = input('Enter the layer for projection: ')
    try:
        layer_index = int(user_input)
        if layer_index > nlayers or layer_index < min_layers:
            sys.exit('FATAL ERROR: Input is out of range. Must be greater '
                     f'than {min_layers-1} and less than {nlayers+1}.')
    except ValueError:
        sys.exit('FATAL ERROR: Input is not a number.')
    return layer_index


def make_energy_array(emin, emax, de, efermi):
    """
    Create the energy array for the plot
    :param emin: minimum energy, ie: ymin (float)
    :param emax: maximum energy, ie: ymax (float)
    :param de: Energy step (float)
    :param efermi: Fermi level (float)
    :return: array of energies: ndarray(dtype=float)
    """
    nene = int((emax - emin) / de)
    energy_window = np.empty(nene)
    for i in range(nene):
        energy_window[i] = (emin - efermi) + (i*de)
    return energy_window


def read_master_input():
    """
    Read the master input file
    :return: nlayer (int), energy_array (1darray(dtype=float)),
             seedname (string), fig_dimensions (1darray(dtype=float))
    """
    nlayers = 0
    plotmode = 1
    fermi_level = 0.0
    emin = 0.0
    emax = 0.0
    energy_step = 1.0
    seedname = 'seedname'
    fig_dimensions = np.empty([2])
    with open('INPUT', 'r', encoding='utf-8') as f:
        for line in f:
            split_line = line.split()
            if split_line[0] == 'nlayers':
                nlayers = int(split_line[1])
            elif split_line[0] == 'e_fermi':
                fermi_level = float(split_line[1])
            elif split_line[0] == 'energy_step':
                energy_step = float(split_line[1])
            elif split_line[0] == 'energy_range':
                emin = float(split_line[1])
                emax = float(split_line[2])
            elif split_line[0] == 'seedname':
                seedname = split_line[1]
            elif split_line[0] == 'figsize':
                fig_dimensions = (float(split_line[1]), float(split_line[2]))
            elif split_line[0] == 'bands_plot':
                band_switch = int(split_line[1])
            elif split_line[0] == 'spectra_plot':
                gf_switch = int(split_line[1])
    if band_switch == 1 and gf_switch == 1:
        plotmode = 3
    elif band_switch == 1 and gf_switch != 1:
        plotmode = 1
    elif band_switch != 1 and gf_switch == 1:
        plotmode = 2
    energy_array = make_energy_array(emin, emax, energy_step, fermi_level)
    return nlayers, energy_array, seedname, fig_dimensions, plotmode, fermi_level


def read_hr(seedname):
    with open(seedname+'_hr.dat', 'r', encoding='utf-8') as f:
        f.readline()
        num_bands = int(f.readline())
    print(f'There are {num_bands} bands.')
    return num_bands


def read_nnkp(seedname):
    """
    Read the nnkp file to get the reciprocal lattice vectors
    :param seedname: the name of the system (string)
    :return: bvec (2darray(dtype=float)) the reciprocal lattice vectors
    """
    filename = seedname + '.nnkp'
    bvec = np.empty([3, 3])
    with open(filename, encoding='utf-8') as f:
        for line in f:
            if line == 'begin recip_lattice\n':
                for i in range(3):
                    bvec[i, :] = f.readline().split()
    return bvec


def read_kpoints(seedname):
    """
    Make the x-axis for the plot.
    :param seedname: the name of the system (string)
    :return: kdists the x-axis for the plot (1darray(dtype=float))
    """
    bvec = read_nnkp(seedname)
    f = open('kpoints', 'r')
    kpt_per_path = int(f.readline())
    n_hsym_pts = int(f.readline())
    hsym_pts = np.empty([n_hsym_pts, 3])
    abs_hsym_pts = np.empty([n_hsym_pts, 3])
    kdists = np.zeros([1 + (kpt_per_path * (n_hsym_pts - 1))])
    for i in range(n_hsym_pts):
        hsym_pts[i, :] = f.readline().split()[0:3]
    f.close()
    for i in range(n_hsym_pts):
        abs_hsym_pts[i, :] = hsym_pts[i, 0] * bvec[0, :] \
                             + hsym_pts[i, 1] * bvec[1, :] \
                             + hsym_pts[i, 2] * bvec[2, :]
    jk=0
    for i in range(n_hsym_pts - 1):
        abs_kpath = abs_hsym_pts[i + 1, :] - abs_hsym_pts[i, :]
        abs_dk = abs_kpath / kpt_per_path
        for ik in range(kpt_per_path):
            kdists[jk + 1] = np.linalg.norm(abs_dk) + kdists[jk]
            jk += 1
    kdists -= kdists[kpt_per_path]
    return kdists


def read_spectral_function(seedname, nene, nkp, nlayers):
    """
    Read the input from the spectral function text file
    :param seedname: the name of the system (string)
    :param nene: number of energies in the y axis
    :param nkp:  number of kpoints in the x axis
    :param nlayers: number of unit-cells used to model the system
    :return: spectral_function (ndarray(dtype=float))
    """
    filename = seedname + '_spec_func.dat'
    spectral_function_1d = np.loadtxt(filename)
    spectral_function = np.reshape(spectral_function_1d, (nlayers, nene, nkp),
                                   order='F')
    return spectral_function


def read_eigenvalues(seedname, nkp):
    """
    Read the plaintext file with the eigenvalue data
    :param seedname: the name of the system (string)
    :param nkp: number of k points - x axis
    :return: eigenvals: (ndarray(dtype=float))
    """
    eigenvals_1d = np.loadtxt(seedname+'_eigenval.dat')
    num_bands = int(len(eigenvals_1d) / nkp)
    eigenvals = np.reshape(eigenvals_1d, (num_bands, nkp), order='F')
    np.savetxt('WS2_eigenval_test', eigenvals)
    return eigenvals


def plot_spectra(layer1, layer2, spectral_function, klist, omegas, seedname,
                 fig_dims):
    """
    Draw the graph and save as raster image.
    :param layer1: First unit cell
    :param layer2: Last unit cell
    :param spectral_function: spectral function - will be plotted
    :param klist: x-axis
    :param omegas: y-axis
    :param seedname: the name of the system (string)
    :param fig_dims: dimensions of the figure
    :return: none
    """
    if layer1 != layer2:
        fig_title = f'Layer = {layer1}-{layer2}'
        filename = f'{seedname}_layer_{layer1}-{layer2}.png'
        plot_func = np.sum(spectral_function[layer1-1:layer2, :, :], axis=0)
    else:
        layer = layer2
        fig_title = f'Layer = {layer}'
        filename = f'{seedname}_layer_{layer}.png'
        plot_func = spectral_function[layer - 1, :, :]
    fig = plt.figure(figsize=fig_dims)
    ax = fig.add_subplot(1, 1, 1)
    yy, xx = np.meshgrid(omegas, klist)
    print('Drawing spectral plot...')
    ax.pcolormesh(xx, yy, plot_func[:, :].T, shading='gouraud',
                  norm='log', cmap='Purples')
    ax.set_title(fig_title)
    ax.set_ylabel(r'$E - E_F$ (eV)')
    ax.set_xlabel(r'$k (\AA^{-1})$')
    plt.tight_layout()
    plt.savefig(filename, dpi=350)


def plot_bands(bandstructure, klist, fig_dims, fermi_level, seedname):
    """
    Plot the band structure
    :param bandstructure: 2d array. Bands are the rows.
    :param klist: the x-axis
    :param fig_dims: dimensions of the figure
    :param fermi_level: the fermi level
    :param seedname: nane of the system
    :return:
    """
    print('Drawing bandstructure plot...')
    fig = plt.figure(figsize=fig_dims)
    ax = fig.add_subplot(1, 1, 1)
    num_bands = bandstructure.shape[0]
    for i in range(num_bands):
        ax.plot(klist, bandstructure[i, :]-fermi_level, color='tab:blue')
    ax.set_ylabel(r'$E - E_F$ (eV)')
    ax.set_xlabel(r'$k (\AA^{-1})$')
    plt.tight_layout()
    plt.savefig(seedname+'_eigenval.pdf')


def main():
    # Initialise all the variables from the input files
    nlayers, energy_array, seedname, fig_dims, plotmode, fermi_level = read_master_input()
    if plotmode == 0:
        sys.exit('You want to plot neither bands nor spectral function...'
                'Quitting the programme...')
    kdists = read_kpoints(seedname)
    nkp = len(kdists)
    nene = len(energy_array)

    if plotmode != 1:
        # Figure out which layers we are going to make a plot for
        layer_index1 = read_layer_index(nlayers, 1)
        layer_index2 = read_layer_index(nlayers, layer_index1)

        print('Reading spectral data...')
        spectral_function = read_spectral_function(seedname, nene, nkp, nlayers)
        plot_spectra(layer_index1, layer_index2, spectral_function, kdists,
                     energy_array, seedname, fig_dims)
    if plotmode != 2:
        print('Reading eigenvalues...')
        eigenvals =  read_eigenvalues(seedname, nkp)
        plot_bands(eigenvals, kdists, fig_dims, fermi_level, seedname)
        


if __name__ == '__main__':
    main()
