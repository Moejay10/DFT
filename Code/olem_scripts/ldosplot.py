#!/usr/bin/env python
# 
from __future__ import division
from __future__ import print_function

from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import os, sys
from subprocess import call

if not os.path.isfile('DOS0'):
    call('split_dos')

if len(sys.argv) < 2:
    print('Usage: ldosplot.py [atom #(s)]')
    exit()

atoms=sys.argv[1:]

for atom in atoms:
    dosfile='DOS'+atom
    DOSdata = loadtxt(dosfile)
    if len(DOSdata[0])==10:
        plot(DOSdata[:,0],DOSdata[:,1],label='s')
        plot(DOSdata[:,0],DOSdata[:,2],label='p')
        plot(DOSdata[:,0],DOSdata[:,5],label='d')
    elif len(DOSdata[0])==4:
        plot(DOSdata[:,0],DOSdata[:,1],label='s')
        plot(DOSdata[:,0],DOSdata[:,2],label='p')
        plot(DOSdata[:,0],DOSdata[:,3],label='d')
    xlabel('$E - E_{Fermi}$ (eV)')
    ylabel('DOS')
    legend()
    suptitle('LDOS atom no '+atom)
    savefig('LDOS'+atom+'.png')

    show()
