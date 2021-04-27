# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 10:35:20 2019

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



# The efficiency vs the thickness.

if __name__ == '__main__':
########### band-gap array #############
    # thickness nm
    #tkArr=np.array([1,2,3,4,5,6,7,8,9,10]) #small
    #tkArr = np.array([100,200,500,1000,2000])
    tkArr=np.array([100,250,500,1000,2000,5000,10000,50000,100000,500000,1000000,10000000,100000000]) # Here we input the thcikness for the code.
        
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

    effDict={} # This is a dictionary.

    for j in range(len(bgArrReal)):
        print('------------------------------------------------------')
        print('------------------------------------------------------')
        effTkArr=[] # A list containing the efficiency data. 
        for i in range(len(tkArr)):
            print("Thickness = %s nm; Band-gap = %s eV; Material = %s " % (str(tkArr[i]),str(bgArrReal[j]),fileName[j]))
            maxEff=calEff(tkArr[i],bgArrReal[j],bgShift[j],0,0,1,fileName[j],bgArrFun[j])
            print(maxEff)
            effTkArr.append(maxEff) 
        effDict[fileName[j]]=effTkArr 

# This makes the ouput file MAXEFF.txt for your result of the efficiency vs the thickness of the material you tested.
File=pd.DataFrame(effDict,index=tkArr)
File.to_csv('MAXEFF.txt', header=True, index=True, sep=str(u'\t'), mode='w')