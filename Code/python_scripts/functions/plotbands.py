#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

import re, os, shutil, sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from ase import io
from ase.dft.kpoints import ibz_points, get_bandpath, special_paths

def band():
    
    # Read crystal structure from POSCAR
    atoms = io.read('POSCAR')

    # Read Fermi level from OUTCAR
    Fermi =   float( os.popen('grep fermi OUTCAR').readlines()[1].split()[2] )

    # Hard code Fermi level from SCF calculation
    Fermi = -0.7460

    # Read band energies from EIGENVAL
    with open('EIGENVAL') as f:
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()
        line4 = f.readline()
        comment = f.readline()
        unknown, npoints, nbands = [int(x) for x in f.readline().split()]
        blankline = f.readline()
        band_energies = [[] for i in range(nbands)]
        for i in range(npoints):
            x, y, z, weight = [float(x) for x in f.readline().split()]
            for j in range(nbands):
                fields = f.readline().split()
                id, energy = int(fields[0]), float(fields[1])
                band_energies[id-1].append(energy)
            blankline = f.readline()
    f.close()
    # Convert to numpy array
    band_energies = np.array(band_energies)
    # Renormalize energy to the Fermi level
    band_energies = band_energies-Fermi

    # Define special points
    # See https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html
    # for special points in different Bravais cells
    points = ibz_points['orthorhombic']
    X = points['X']
    G = points['Gamma']
    Y = points['Y']
    T = points['T']
    Z = points['Z']


    # Define path to be plotted. Has to be the same that was 
    # calculated by VASP (defined in KPOINTS)
    path = get_bandpath([G, Y, T, Z, G], atoms.cell, npoints=npoints)
    x2, X2, labels = path.get_linear_kpoint_axis() 
    labels = [label.replace('G','\Gamma') for label in labels]


    # Hardcode plots with even distance between special points:
    nkpoints=int(npoints/(len(labels)-1))
    x1 = list(range(npoints))
    X1 = [x1[0]]+x1[nkpoints-1::nkpoints]

    xticks(X1, ['$%s$' % n for n in labels])

    # Create plot
    for iband in range(nbands):
        plot(x1,band_energies[iband,:],color='k')


    # Change range:
        ylim(-2,2)
    
    # Making the line for zero point
    a = linspace(1,100)
    b = zeros(len(a))
    plot(a,b,"--k")


    xlabel('k points')
    ylabel('$E - E_{Fermi}$ [eV]')

    # Print figure to file
    #savefig('BaSi2.png')

    show()

if __name__='__main__':
    band()