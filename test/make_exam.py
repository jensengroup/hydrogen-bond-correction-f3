#!/usr/bin/env python

import numpy as np
import subprocess
import csv
import re

# S66
s66_hydrogenb  = range(1, 24)  # 1 .. 23
s66_dispersion = range(24, 47) # 24 .. 46

######################################################################

def return_code(execute):
    """ Return the first float of a bash command
    """
    string = subprocess.Popen(execute, shell=True, stdout=subprocess.PIPE).communicate()[0]
    number = re.findall(r'[-]*\d+\.\d+', string)
    number = number[0]
    number = float(number)
    return number


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

if __name__ == '__main__':

    s66 = []
    s66_lib = []
    s66_pm6 = []
    s66_hb = []

    s66_model = []

    # S66
    cr = csv.reader(open('begdb_s66_energies.csv', 'rb'), quotechar="'")

    for row in cr:
        sys_id = row[0]
        sys_name = row[1]
        sys_energy = row[3]
        sys_energy = float(sys_energy)

        sys_filename = sys_id+sys_name
        sys_filename = sys_filename.replace('...', '')
        sys_filename = sys_filename.replace('(', '')
        sys_filename = sys_filename.replace(')', '')
        sys_filename = sys_filename.replace(' ', '')
        sys_filename = sys_filename.replace('-', '')

        s66_lib.append(sys_energy)
        s66.append(sys_filename)


    # PM6-D3
    cr = csv.reader(open('pm6-d3_energies.csv', 'rb'))

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

    for x in s66:
        energy_a = return_code('../f3_exe -p param.dat structures/'+x+'_a.xyz | grep -A3 "Correction Energy" | grep kcal')
        energy_b = return_code('../f3_exe -p param.dat structures/'+x+'_b.xyz | grep -A3 "Correction Energy" | grep kcal')
        energy_c = return_code('../f3_exe -p param.dat structures/'+x+'_c.xyz | grep -A3 "Correction Energy" | grep kcal')
        s66_hb.append(energy_c - energy_a - energy_b)


    # FULL RMSD
    for i in range(len(s66)):
        # print "{0:26s} {1:10.4f} {2:10.4f} {3:10.4f} ".format(s66[i], s66_lib[i], s66_pm6[i]+s66_hb[i], s66_pm6[i]+s66_hb[i]-s66_lib[i])
        s66_model.append(s66_pm6[i]+s66_hb[i])

    s66_model = np.array(s66_model)
    s66_lib = np.array(s66_lib)

    s66_rmsd = rmsd(s66_model, s66_lib)

    # print "RMSD: {0:10.4f}".format(s66_rmsd)
    # print "MAX:  {0:10.4f}".format(max(abs(s66_lib-s66_model)))


    # HYDROGEN BOND RMSD
    s66_model = []

    for i in range(len(s66_hydrogenb)):
        print "{0:26s} {1:10.4f} {2:10.4f} {3:10.4f} ".format(s66[i], s66_lib[i], s66_pm6[i]+s66_hb[i], s66_pm6[i]+s66_hb[i]-s66_lib[i])
        s66_model.append(s66_pm6[i]+s66_hb[i])

    s66_model = np.array(s66_model)
    s66_rmsd = rmsd(s66_model, s66_lib[:23])
    print "RMSD: {0:10.4f}".format(s66_rmsd)
    print "MAX:  {0:10.4f}".format(max(abs(s66_lib[:23]-s66_model)))


