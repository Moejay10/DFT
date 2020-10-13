#!/usr/bin/env python

from ase.io import read, vasp
from ase.build import sort, surface

unitcell = sort(read('POSCAR'))

slab = surface(unitcell, indices=(1,1,0), layers=4, vacuum=10)

vasp.write_vasp('POSCAR_slab',sort(slab),vasp5=True,direct=True)
