#!/usr/bin/env python

import numpy as np
import subprocess
import csv
import re

import multiprocessing as mp

#MATPLOTLIB
import matplotlib
matplotlib.use('Agg')

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from matplotlib.image import NonUniformImage
import matplotlib.colors as col

# Use Latex
from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('xtick', labelsize=15)
rc('ytick', labelsize=15)
matplotlib.rcParams.update({'font.size': 20})


# S66
s66_hb   = range(1, 24)  # 1 .. 23
s66_disp = range(24, 47) # 24 .. 46

# S22
s22_hb   = range(1, 8)  # 1 .. 7
s22_disp = range(8, 16) # 8 .. 15


# How many CPU's can I use?
WORKERS = mp.cpu_count()

######################################################################

def return_code(execute):
    """ Return the first float of a bash command
    """
    string = subprocess.Popen(execute, shell=True, stdout=subprocess.PIPE).communicate()[0]
    number = re.findall(r'[-]*\d+\.\d+', string)
    number = number[0]
    number = float(number)
    return number


def split_jobs(jobs, workers):
  """ Split job array into MP
  """
  return [jobs[i::workers] for i in xrange(workers)]


def mad(L1, L2):
    """ Mean Absolute Deviation
        MAD = 1/n * \sum_{i=1}^n abs( x_i - \hat x )
    """
    LD = np.abs(L1 - L2)
    mad = sum(LD)/len(LD)
    return mad


def md(L1, L2):
    """ Mean deviation
    """
    LD = L1 - L2
    md = sum(LD)/len(LD)
    return md


def rmsd(L1, L2):
    """ Caculate Root-mean-square deviation (RMSD) from two vectors
    """
    LD = (L1 - L2)**2
    lrmsd = sum(LD)
    lrmsd = np.sqrt(lrmsd/len(LD))
    return lrmsd


######################################################################


def calculate_hydrogen(param, structure_set):
    """ Calculate hydrogen correction given a set of parameters, param
    """
    corr_energies = []

    c_nitrogen, c_oxygen = param
    filename = str(c_nitrogen)+'_'+str(c_oxygen)
    folder = 'parameters'

    # Create param file
    f = open(folder+'/'+filename, 'w')
    f.write(str(c_nitrogen))
    f.write('\n')
    f.write(str(c_oxygen))
    f.close()

    for x in structure_set:
        energy_a = return_code('../f3_exe -p '+folder+'/'+filename+' structures/'+x+'_a.xyz | grep -A4 "Correction" | grep kcal')
        energy_b = return_code('../f3_exe -p '+folder+'/'+filename+' structures/'+x+'_b.xyz | grep -A4 "Correction" | grep kcal')
        energy_c = return_code('../f3_exe -p '+folder+'/'+filename+' structures/'+x+'_c.xyz | grep -A4 "Correction" | grep kcal')

        corr_energies.append(energy_c - energy_a - energy_b)

    return np.array(corr_energies)

######################################################################

if __name__ == '__main__':

    s22 = []
    s66 = []
    s22_lib = []
    s66_lib = []
    s22_pm6 = []
    s66_pm6 = []

    # S22
    cr = csv.reader(open('energies_s22.csv', 'rb'), quotechar="'")
    for row in cr:
        sys_id = row[0]
        sys_name = row[1]
        sys_energy = row[3]
        sys_energy = float(sys_energy)

        sys_filename = 's22/'+sys_id+sys_name
        sys_filename = sys_filename.replace('...', '')
        sys_filename = sys_filename.replace('(', '')
        sys_filename = sys_filename.replace(')', '')
        sys_filename = sys_filename.replace(' ', '')
        sys_filename = sys_filename.replace('-', '')

        s22_lib.append(sys_energy)
        s22.append(sys_filename)


    # S66
    cr = csv.reader(open('energies_s66.csv', 'rb'), quotechar="'")
    for row in cr:
        sys_id = row[0]
        sys_name = row[1]
        sys_energy = row[3]
        sys_energy = float(sys_energy)

        sys_filename = 's66/'+sys_id+sys_name
        sys_filename = sys_filename.replace('...', '')
        sys_filename = sys_filename.replace('(', '')
        sys_filename = sys_filename.replace(')', '')
        sys_filename = sys_filename.replace(' ', '')
        sys_filename = sys_filename.replace('-', '')

        s66_lib.append(sys_energy)
        s66.append(sys_filename)


    # S22 PM6-D3
    cr = csv.reader(open('energies_s22_pm6-d3.csv', 'rb'))
    for row in cr:
        sys_id = row[0]
        sys_name = row[1]
        sys_energy = row[2]
        sys_energy = float(sys_energy)

        sys_filename = sys_id+sys_name
        sys_filename = sys_filename.replace('...', '')
        sys_filename = sys_filename.replace('(', '')
        sys_filename = sys_filename.replace(')', '')
        sys_filename = sys_filename.replace(' ', '')
        sys_filename = sys_filename.replace('-', '')

        s22_pm6.append(sys_energy)


    # S66 PM6-D3
    cr = csv.reader(open('energies_s66_pm6-d3.csv', 'rb'))
    for row in cr:
        sys_id = row[0]
        sys_name = row[1]
        sys_energy = row[2]
        sys_energy = float(sys_energy)

        sys_filename = sys_id+sys_name
        sys_filename = sys_filename.replace('...', '')
        sys_filename = sys_filename.replace('(', '')
        sys_filename = sys_filename.replace(')', '')
        sys_filename = sys_filename.replace(' ', '')
        sys_filename = sys_filename.replace('-', '')

        s66_pm6.append(sys_energy)


    shb = s22[:7] + s66[:23]
    shb_lib = s22_lib[:7] + s66_lib[:23]
    shb_pm6 = s22_pm6[:7] + s66_pm6[:23]


    # SCAN PARAMETER SPACE
    # should be around ~-0.1
    scan_nitrogen = np.arange(0.0, -0.2, -0.01)
    scan_oxygen   = np.arange(0.0, -0.2, -0.01)

    n_nitrogen = len(scan_nitrogen)
    n_oxygen   = len(scan_oxygen)
    n_total = n_nitrogen*n_oxygen

    full_rmsd = np.zeros((n_nitrogen, n_oxygen))
    S = full_rmsd

    print "Scanning", n_total, " energies", "using", WORKERS, 'cpu\'s'

    # for i in xrange(n_nitrogen):
    #     for j in xrange(n_oxygen):

    # TODO
    # http://briansimulator.org/sharing-numpy-arrays-between-processes
    from multiprocessing import sharedctypes
    size = S.size
    shape = S.shape
    S.shape = size
    S_ctypes = sharedctypes.RawArray('d', S)
    S = np.frombuffer(S_ctypes, dtype=np.float64, count=size)
    S.shape = shape

    from numpy import ctypeslib

    def worker(id, job):
        """ worker function for MP """

        S = ctypeslib.as_array(S_ctypes)
        S.shape = shape

        for i in job:
            for j in xrange(n_oxygen):

                N = scan_nitrogen[i]
                O = scan_oxygen[j]
                param = (N, O)

                hb_energies = calculate_hydrogen(param, shb)
                param_energies = shb_pm6 + hb_energies
                param_rmsd = rmsd(shb_lib, param_energies)
                S[i, j] = param_rmsd


    job_list = split_jobs(range(n_nitrogen), WORKERS)

    processes = [mp.Process(target=worker, args=(i, job_list[i] )) for i in xrange(WORKERS)]
    for p in processes: p.start() # Fork
    for p in processes: p.join()  # Join

    full_rmsd = S

    flat_idx = full_rmsd.argmin()
    i, j = np.unravel_index(flat_idx, full_rmsd.shape)

    print 'c_nitrogen: ', scan_nitrogen[i]
    print 'c_oxygen:   ', scan_oxygen[j]
    print 'rmsd min    ', np.min(full_rmsd)

    extent = [scan_nitrogen[0], scan_nitrogen[-1], scan_oxygen[0], scan_oxygen[-1]]

    im = plt.imshow(full_rmsd,
                    interpolation='nearest',
                    extent=extent,
                    cmap='Spectral')

    # Labels
    plt.xlabel('$\mathrm{C}_\mathrm{O}$')
    plt.ylabel('$\mathrm{C}_\mathrm{N}$')
    # plt.title('Scan of NSP and OSP')

    # Color-bar
    cbar = plt.colorbar(im)
    cbar.set_label('$\mathrm{RMSD [kcal/mol]}$', rotation=270)

    # Grid
    plt.grid(True)

    # Save figure
    plt.savefig('hydrogen_rmsd_scan.png')
    plt.savefig('hydrogen_rmsd_scan.tif')
    plt.savefig('hydrogen_rmsd_scan.eps')

