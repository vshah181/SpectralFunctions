import numpy as np
import sys
import matplotlib.pyplot as plt


def read_layer_index(nlayers, min_layers):
    user_input = input('Enter the layer for projection: ')
                       
    try:
        layer_index = int(user_input)
        if layer_index > nlayers or layer_index < min_layers:
            sys.exit('FATAL ERROR: Input is out of range. Must be greater '\
                    f'than {min_layers-1} and less than {nlayers+1}.')
    except ValueError:
        sys.exit('FATAL ERROR: Input is not a number.')
    return layer_index


def make_energy_array(emin, emax, de, efermi):
    nene = int((emax - emin) / de)
    energy_window = np.empty(nene)
    for i in range(nene):
        energy_window[i] = (emin - efermi) + (i*de)
    return energy_window


def read_master_input():
    f = open('INPUT', 'r')
    nlayers = 0
    fermi_level = 0.0
    emin = 0.0
    emax = 0.0
    energy_step = 1.0
    seedname = 'seedname'
    fig_dimensions = np.empty([2])
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
    f.close()
    energy_array = make_energy_array(emin, emax, energy_step, fermi_level)
    return nlayers, fermi_level, energy_array, seedname, fig_dimensions


def read_nnkp(seedname):
    filename = seedname + '.nnkp'
    bvec = np.empty([3, 3])
    f = open(filename)
    for line in f:
        if line == 'begin recip_lattice\n':
            for i in range(3):
                bvec[i, :] = f.readline().split()
    f.close()
    return bvec


def read_kpoints(seedname):
    f = open('kpoints', 'r')
    bvec = read_nnkp(seedname)
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
    ik = 0
    for i in range(n_hsym_pts - 1):
        kpath = hsym_pts[i + 1, :] - hsym_pts[i, :]
        dk = kpath / kpt_per_path
        abs_kpath = abs_hsym_pts[i + 1, :] - abs_hsym_pts[i, :]
        abs_dk = abs_kpath / kpt_per_path
        for j in range(kpt_per_path):
            kdists[ik + 1] = np.linalg.norm(abs_dk) + kdists[ik]
            ik += 1
    kdists -= kdists[kpt_per_path]
    return kdists


def read_spectral_function(seedname, nene, nkp, nlayers):
    filename = seedname + '_spec_func.dat'
    spectral_function_1d = np.loadtxt(filename)
    spectral_function = np.reshape(spectral_function_1d, (nlayers, nene, nkp),
                                   order='F')
    return spectral_function


def plot_spectra(layer1, layer2, spectral_function, klist, omegas, seedname,
                 fig_dims):
    if layer1 != layer2:
        fig_title = f'Layer = {layer1}-{layer2}'
        filename = f'{seedname}_layer_{layer1}-{layer2}.png'
        plot_func = (np.sum(spectral_function[layer1-1:layer2, :, :], axis=0))
    else:
        layer = layer2
        fig_title = f'Layer = {layer}'
        filename = f'{seedname}_layer_{layer}.png'
        plot_func = spectral_function[layer - 1, :, :]
    fig = plt.figure(figsize=fig_dims)
    ax = fig.add_subplot(1, 1, 1)
    yy, xx = np.meshgrid(omegas, klist)
    print('Drawing plot...')
    ax.pcolormesh(xx, yy, plot_func[:, :].T, shading='gouraud',
                  norm='log', cmap='Purples')
    ax.set_title(fig_title)
    ax.set_ylabel(r'$E - E_F$ (eV)')
    ax.set_xlabel(r'$k (\AA^{-1})$')
    plt.tight_layout()
    plt.savefig(filename, dpi=350)


def main():
    # Initialise all the variables from the input files
    nlayers, efermi, energy_array, seedname, fig_dims = read_master_input()
    kdists = read_kpoints(seedname)
    nkp = len(kdists)
    nene = len(energy_array)

    print('Reading input file...')
    spectral_function = read_spectral_function(seedname, nene, nkp, nlayers)

    # Figure out which layers we are going to make a plot for
    layer_index1 = read_layer_index(nlayers, 1)
    layer_index2 = read_layer_index(nlayers, layer_index1)
    plot_spectra(layer_index1, layer_index2, spectral_function, kdists,
                 energy_array, seedname, fig_dims)


if __name__ == '__main__':
    main()
