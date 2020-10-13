#!/usr/bin/env python

from ase.io import read,vasp
from ase.build import sort, cut

unitcell = sort(read('POSCAR'))

supercell = cut(unitcell,a=(2,0,0),b=(0,2,0),c=(0,0,2))

vasp.write_vasp('POSCAR_supercell',sort(supercell),vasp5=True,direct=True)
