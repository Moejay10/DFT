#!/bin/env python
#
# Fit output data from VASP to the Birch-Murnaghan
# equation of state.
# Written by Ole Martin Lovvik, 2018-09-17
# $LastChangedDate$
# $Rev$
#
# Please send bugs to ole.martin.lovvik@sintef.no
'''
Usage: fit_birch_murnaghan.py [options] [inputfile]

An inputfile has to be present in the starting folder.
If not specified, the script will look for "results.dat".
The file should be generated with the command
toten */OUTCAR > results.dat
and have the following format:

Comment_line
a1 F1 E01
a2 F2 E02
...

a1, a2, ...   : lattice constants
F1, F2, ...   : free energies
E01, E02, ... : total energies

The fitted curve is plotted to the file "birch_murnaghan.png"

Options:
  -h, --help              Show this help message and exit
'''

from __future__ import division
from __future__ import print_function

import sys, getopt
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

shopts = 'h'
longopts = ['help']

try:
   opts, args = getopt.getopt(sys.argv[1:], shopts, longopts)
except getopt.GetoptError as err:
   # print help information and exit:
   print('{0}: {1}'.format(sys.argv[0], str(err)))
   print(__doc__)
   sys.exit(2)
   
for o, a in opts:
   if o in ('-h', '--help'):
      # print help information and exit:
      print(__doc__)
      sys.exit()
   else:
      assert False, 'unhandled option'
      print('Error: option is not known')
      sys.exit(2)

if len(args) == 0:
   inputfile = "results.dat"
elif len(args) == 1:
   inputfile = args[0]
else:
   # print help information and exit:
   print('{0}: {1}'.format(sys.argv[0], str(err)))
   print(__doc__)
   sys.exit(2)
   
# Birch-Murnaghan equation of state (https://en.wikipedia.org/wiki/Birch%E2%80%93Murnaghan_equation_of_state):
def Birch_Murnaghan(a,V0,E0,B0,B0p):
   # a: lattice constant
   # V0: Equilibrium volume
   # E0: Equilibrium energy
   # B0: Bulk modulus
   # B0p: p-derivative of B0
   # Assume cubic unit cell:
   V = a**3
   return ( E0 + (9*V0*B0)/16 * (\
         ((V0/V)**(2/3) - 1)**3 * B0p +\
            ((V0/V)**(2/3) - 1)**2 *\
            (6 - 4*(V0/V)**2/3)))
                  
# Read data:
try:
   adata,efree,toten = np.loadtxt(inputfile,unpack=True, skiprows=1)
except:
   print("Trouble reading data from inputfile")
   print(__doc__)
   exit(1)

# Select middle value 
def findMiddle(input_list):
    middle = float(len(input_list))/2
    if middle % 2 != 0:
        return input_list[int(middle - .5)]
    else:
        return (input_list[int(middle)], input_list[int(middle-1)])

a_init=findMiddle(adata)
e_init=findMiddle(toten)

# Initial guess:
a0=a_init    #Å 
V00=a0**3 #Å**3
E00=e_init   #eV
B00=0.5    #eV/Å**3
B0p0=1    #1
p0=[V00,E00,B00,B0p0]

# Units:
eV=1.6021766208e-19 # J
AA=1e-10 # m
GPa=eV/(AA**3)*1e-9 # Pa

popt,pcov = curve_fit(Birch_Murnaghan,adata,toten,p0)

#print("Optimized parameters: ",popt)
p_sigma = np.sqrt(np.diag(pcov))
#print(p_sigma)

print("Results from the fitting:")
print("--------------------------------------------")
print("Optimized volume (AA^3):          ",popt[0])
print("Optimized lattice constant (AA):  ",float(popt[0])**(1/3))
print("Optimized energy (eV):            ",popt[1])
print("Optimized bulk modulus (eV/AA^3): ",popt[2])
print("Optimized bulk modulus (GPa):     ",float(popt[2])*GPa)

#exit()
a_plot = np.linspace(adata[0],adata[-1],50)
#toten_fit = popt[0]+popt[1]*(a_plot-popt[2])**2


plt.plot(a_plot, Birch_Murnaghan(a_plot, *popt), 'k-')
plt.plot(adata,toten,'o')
plt.xlabel('Lattice const. (AA)')#, fontsize=14)
plt.ylabel('Total energy (eV)')
plt.savefig("birch_murnaghan.png")
plt.show()
