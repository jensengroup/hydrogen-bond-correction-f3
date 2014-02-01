#!/usr/bin/env python

import numpy as np
import subprocess
import csv
import re

# S66
s66_hb   = range(1, 24)  # 1 .. 23
s66_disp = range(24, 47) # 24 .. 46

# S22
s22_hb   = range(1, 8)  # 1 .. 7
s22_disp = range(8, 16) # 8 .. 15

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
    s22 = []

    s66_lib = []
    s22_lib = []

    s66_model = []
    s22_model = []

    # S22
    cr = csv.reader(open('energies_s22.csv', 'rb'), quotechar="'")
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

        s22_lib.append(sys_energy)
        s22.append(sys_filename)


    # S66
    cr = csv.reader(open('energies_s66.csv', 'rb'), quotechar="'")
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


    # S22 PM6-D3
    cr = csv.reader(open('energies_s22_pm6-d3_2.csv', 'rb'))
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

        x = sys_filename

        energy_a = return_code('../f3_exe -p param.dat structures/s22/'+x+'_a.xyz | grep -A3 "Correction Energy" | grep kcal')
        energy_b = return_code('../f3_exe -p param.dat structures/s22/'+x+'_b.xyz | grep -A3 "Correction Energy" | grep kcal')
        energy_c = return_code('../f3_exe -p param.dat structures/s22/'+x+'_c.xyz | grep -A3 "Correction Energy" | grep kcal')

        # PM6-D3H+
        s22_model.append(sys_energy + energy_c - energy_a - energy_b)


    # S66 PM6-D3
    cr = csv.reader(open('energies_s66_pm6-d3_2.csv', 'rb'))
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

        x = sys_filename

        energy_a = return_code('../f3_exe -p param.dat structures/s66/'+x+'_a.xyz | grep -A3 "Correction Energy" | grep kcal')
        energy_b = return_code('../f3_exe -p param.dat structures/s66/'+x+'_b.xyz | grep -A3 "Correction Energy" | grep kcal')
        energy_c = return_code('../f3_exe -p param.dat structures/s66/'+x+'_c.xyz | grep -A3 "Correction Energy" | grep kcal')

        # PM6-D3H+
        s66_model.append(sys_energy + energy_c - energy_a - energy_b)


    sxx_model = np.array(s22_model + s66_model)
    sxx_lib   = np.array(s22_lib   + s66_lib)
    sxx_rmsd = rmsd(sxx_model, sxx_lib)

    shb_model = []
    shb_lib = []

    print "Full sets"
    print "RMSD: {0:10.4f}".format(sxx_rmsd)
    print "MAX:  {0:10.4f}".format(max(abs(sxx_lib-sxx_model)))
    print
    print "Hydrogen-Bonds"
    print
    print "S22"

    for i in s22_hb:
        j = i -1
        print " {0:26s} {1:10.4f} {2:10.4f} {3:10.4f} ".format(s22[j], s22_lib[j], s22_model[j], s22_model[j]-s22_lib[j])
        shb_model.append(s22_model[j])
        shb_lib.append(s22_lib[j])

    print
    print "S66"

    for i in s66_hb:
        j = i - 1
        print " {0:26s} {1:10.4f} {2:10.4f} {3:10.4f} ".format(s66[j], s66_lib[j], s66_model[j], s66_model[j]-s66_lib[j])
        shb_model.append(s66_model[j])
        shb_lib.append(s66_lib[j])

    print
    shb_model = np.array(shb_model)
    shb_lib = np.array(shb_lib)
    shb_rmsd = rmsd(shb_model, shb_lib)
    print "RMSD: {0:10.4f}".format(shb_rmsd)
    print "MAX:  {0:10.4f}".format(max(abs(shb_lib-shb_model)))

