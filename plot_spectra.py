import numpy as np
import matplotlib.pyplot as plt


def read_master_input():
    f = open('INPUT', 'r')
    nlayers = 0
    fermi_level = 0.0
    emin = 0.0
    emax = 0.0
    energy_step = 1.0
    case = 'seedname'
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
            case = split_line[1]
    f.close()
    energy_array = np.arange(emin-fermi_level, emax-fermi_level, energy_step)
    return nlayers, fermi_level, energy_array, case


def read_nnkp(case):
    filename = case + '.nnkp'
    bvec = np.empty([3, 3])
    f = open(filename)
    for line in f:
        if line == 'begin recip_lattice\n':
            for i in range(3):
                bvec[i, :] = f.readline().split()
    f.close()
    return bvec


def read_kpoints(case):
    f = open('kpoints', 'r')
    bvec = read_nnkp(case)
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
        kpath = hsym_pts[i+1, :] - hsym_pts[i, :]
        dk = kpath / kpt_per_path
        abs_kpath = abs_hsym_pts[i+1, :] - abs_hsym_pts[i, :]
        abs_dk = abs_kpath / kpt_per_path
        for j in range(kpt_per_path):
            kdists[ik + 1] = np.linalg.norm(abs_dk) + kdists[ik]
            ik += 1
    return kdists


def read_spectral_function(case, nene, nkp, nlayers):
    spectral_function = np.empty([nkp, 1+nlayers, nene])
    filename = case + '_spec_func.dat'
    layer = 1
    f = open(filename, 'r')
    for line in f:
        empty_line = (line == ' \n')
        if not(empty_line) and line.split()[0] == 'layer=':
            try:
                layer = int(line.split()[1])
            except ValueError:
                if line.split()[1] == 'all':
                    layer = 0
            for i in range(nene):
                spectral_function[:, layer, i] = f.readline().split()
    f.close()
    return spectral_function


def plot_spectra(layer, spectral_function, klist, omegas):
    # Layer 0 means plot the sum of all the layer contribution,
    # otherwise you can specify which layer you want to plot.
    fig = plt.figure(figsize=(3, 6))
    ax = fig.add_subplot(1, 1, 1)
    yy, xx = np.meshgrid(omegas, klist)
    ax.pcolormesh(xx, yy, spectral_function[:, layer, :], shading='gouraud',
                  norm='log', cmap='afmhot')
    plt.show()


def main():
    nlayers, efermi, energy_array, case = read_master_input()
    kdists = read_kpoints(case)
    nkp = len(kdists)
    nene = len(energy_array)
    spectral_function = read_spectral_function(case, nene, nkp, nlayers)
    plot_spectra(0, spectral_function, kdists, energy_array)


if __name__ == '__main__':
    main()
