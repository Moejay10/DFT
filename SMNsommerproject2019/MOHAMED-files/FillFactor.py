# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 10:37:01 2019

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


# This code gives the short circuit current density, the open circuit voltage and the fill factor. 

########### band-gap array #############

#This lets you decide what the thickness should be.
# thickness nm
#tkArr=np.array([1,2,3,4,5,6,7,8,9,10]) #small
#tkArr=np.array([5,10,25,50,100,250,500,1000,2000,5000,10000,50000,100000,500000,1000000]) #big
#tkArr = np.array([10,100,1000,10000])
tkArr = np.array([100])


# This lets you decide what the bandgaps should be.
# band gap of material eV
Eg= np.linspace(0.4,4.0,10) # This program is slow, so try not to calculate many points. 
#Eg = np.array([1.0])


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
    fileName = [str(input("fileName="))] 

        
bgArrReal=np.array([DBgMat]*len(fileName)) # minimum direct bandgap
bgArrFun=bgArrReal
bgShift=np.array([IBgMat]*len(fileName)) # minimum bandgap (can be indirect)


Voc = [] # This is a list that contains the voltage data.
Jsc = [] # This is a list that contains the current density data.
FF = [] # This is a list that contains the fillfactor data.

for j in range(len(tkArr)):
    print('------------------------------------------------------')
    print('------------------------------------------------------')
    for i in range(len(Eg)):
        print("Thickness = %s nm; Band-gap = %s eV; Material = %s " % (str(tkArr[j]),str(Eg[i]),fileName[j]))
        maxEff=nvsEg_FF(tkArr[j],bgArrReal[j],bgShift[j],0,0,1,fileName[j],bgArrFun[j],Eg[i])
        print(maxEff[3])
        Voc.append(maxEff[1])
        Jsc.append(maxEff[2])
        FF.append(maxEff[3]) 
    
    



# This makes the ouput file FillFactor.txt for your result of the Voc,Jsc and FF vs the bandgap for different thicknesses of the material you tested.
outfile = open("FillFactor.txt", 'w')
outfile.write("Voc")
outfile.write('\n')
for x in Voc:
    outfile.write('%g' % x)
    outfile.write('\n')
outfile.write('\n')
outfile.write("Jsc")
outfile.write('\n')
for y in Jsc:
    outfile.write('%g' % y)
    outfile.write('\n')
outfile.write('\n')
outfile.write("FF")
outfile.write('\n')
for z in FF:
    outfile.write('%g' % z)
    outfile.write('\n')
outfile.close()