#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

import re, sys
import ase.io.vasp
import ase.io.cif


try:
    infilename  = sys.argv[1]
except:
    # try block failed,
    print("Usage:", sys.argv[0], "infile.cif")
    sys.exit(1)

filename = re.sub(r'.cif', '', infilename)

a = ase.io.read(infilename)
ase.io.vasp.write_vasp('POSCAR',a,label=filename,vasp5=True,direct=True)
