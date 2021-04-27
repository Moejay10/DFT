# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 10:36:25 2019

@author: msismail
"""
from __future__ import division, print_function, unicode_literals
import numpy, scipy.interpolate, scipy.integrate
from pylab import *
import matplotlib.pyplot as plt
import bigfloat
import time

import pandas as pd

import sys

from cal_Eff_FF import *


# The efficiency vs the bandgap

# Inputs for the material, bandgaps (direct & indirect) of material eV and the file with the absorption coefficient data. Must type in this information on the command line when running the code.
try:
    Material = str(sys.argv[1]) # This lets you input on the command line which material you are running the code on.
    DBgMat = float(sys.argv[2]) # This lets you input on the command line  what the direct bandgap of your material is.
    IBgMat = float(sys.argv[3]) # This lets you input on the command line  what the indirect bandap of your material is.
    fileName = [str(sys.argv[4])] # This lets you input on the command line which absorption coefficient file you are using, which then depends on which material you are using.
    
    
except:
    print ("The material, bandgaps (direct and indirect) and the absorption coefficient file must be inputed on the command line") # If you did not type any of the variables on the command line you get a second chance to do so.
    Material = str(input("Material="))
    DBgMat = float(input("Eg_dir="))
    IBgMat = float(input("Eg_ind="))
    fileName = str(input("fileName=")) 
    


bgArrReal= DBgMat # minimum direct bandgap
bgArrFun = bgArrReal
bgShift = IBgMat # minimum bandgap (can be indirect)


# This lets you decide what the bandgaps should be.
# bandgap
Eg= np.linspace(0.2,2.0,10)  # This program is slow, so try not to calculate many points. 


#This lets you decide what the thickness should be
# thickness 
tkArr = np.array([10,100,1000,10000])
#tkArr = np.array([100000000]) # This is a test to see what happens when the thickness is "infinity".


effDict={} # This is a dictionary.

for j in range(len(tkArr)):
        print('------------------------------------------------------')
        print('------------------------------------------------------')
        effTkArr=[] # A list containing the efficiency data.
        for i in range(len(Eg)):
            print("Bandgap = %s eV; Thickness = %s nm; Material = %s " % (str(Eg[i]),str(tkArr[j]),fileName))
            maxEff=nvsEg_FF(tkArr[j],bgArrReal,bgShift,0,0,1,fileName,bgArrFun,Eg[i])[0]
            print(maxEff)
            effTkArr.append(maxEff)
        effDict[fileName[j]]=effTkArr

# This makes the ouput file EFFvsEg.txt for your result of the efficiency vs the bandgap of the material you tested.
File=pd.DataFrame(effDict,index=Eg)
File.to_csv('EFFvsEg.txt', header=True, index=True, sep=str(u'\t'), mode='w')
