# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 00:10:44 2019

@author: mohe9
"""

from __future__ import division, print_function, unicode_literals
import numpy, scipy.interpolate, scipy.integrate
from pylab import *
import matplotlib.pyplot as plt
import bigfloat
import time

import pandas as pd

import sys

from calEffM import *


"""
#Original

if __name__ == '__main__':
########### band-gap array #############
    # thickness nm
    #tkArr=np.array([1,2,3,4,5,6,7,8,9,10]) #small
    tkArr=np.array([1,5,10,25,50,100,250,500,1000,2000,5000]) #big
        
    # band gap of material eV
    #BgCISe=0.848945
    #BgMat = BgCISe
    #BgCIS = 0.87
    #BgMat = BgCIS
    BgGaAs= 1.441
    BgMat = BgGaAs  
    #BgCGS= 2.0
    #BgMat = BgCGS   
    
    # thickness nm
    #tkArr=np.array([30, 100, 500, 2000])
    #tkArr=np.linspace(100, 2000, num=11)

    #fileName=['K12.txt','K16.txt','K18.txt','K20.txt','K24.txt','K30.txt']
    # better smooth out the value and starts from exact bandgap engergy
    #fileName=['CISe_ABS.txt']
    #fileName=['CIS_ABS.txt']
    fileName=['GaAs_ABS.txt']
    #fileName=['CGS_ABS.txt']
        
    bgArrReal=np.array([BgMat]*len(fileName)) # minimum direct bandgap
    bgArrFun=bgArrReal
    bgShift=np.array([BgMat]*len(fileName)) # minimum bandgap (can be indirect)

    effDict={}

    for j in range(len(bgArrReal)):
        print('------------------------------------------------------')
        print('------------------------------------------------------')
        effTkArr=[]
        for i in range(len(tkArr)):
            print("Thickness = %s nm; Band-gap = %s eV; Material = %s " % (str(tkArr[i]),str(bgArrReal[j]),fileName[j]))
            maxEff=calEff(tkArr[i],bgArrReal[j],bgShift[j],0,0,1,fileName[j],bgArrFun[j])
            print(maxEff)
            effTkArr.append(maxEff) 
        effDict[fileName[j]]=effTkArr 
    
outfile = open("test.txt", 'w')
for row in tkArr:
    outfile.write('%0.1f' % row)
    outfile.write('\n')
outfile.write('\n')
for column in effTkArr:
    outfile.write('%.13f' % column)
    outfile.write('\n')
outfile.close()


File=pd.DataFrame(effDict,index=tkArr)
File.to_csv('MAXEFF.txt', header=True, index=True, sep=str(u'\t'), mode='w')
"""



# n vs d 

if __name__ == '__main__':
########### band-gap array #############
    # thickness nm
    #tkArr=np.array([1,2,3,4,5,6,7,8,9,10]) #small
    #tkArr = np.array([100,200,500,1000,2000])
    tkArr=np.array([100,250,500,1000,2000,5000,10000,50000,100000,500000,1000000,10000000,100000000]) #big
        
    # band gap of material eV
    Egap = 0.3 # eV
    try:
      Material = str(sys.argv[1])
      DBgMat = float(sys.argv[2])
      IBgMat = float(sys.argv[3])
      fileName = [str(sys.argv[4])]
    
    except:
      print ("The Material and the bandgaps must be inputed on the command line")
      Material = str(input("Material="))
      DBgMat = float(input("Eg_dir="))
      IBgMat = float(input("Eg_ind=")) 
      fileName = [str(input("fileName="))] 

        
    bgArrReal=np.array([DBgMat]*len(fileName)) # minimum direct bandgap
    bgArrFun=bgArrReal
    bgShift=np.array([IBgMat]*len(fileName)) # minimum bandgap (can be indirect)

    effDict={}

    for j in range(len(bgArrReal)):
        print('------------------------------------------------------')
        print('------------------------------------------------------')
        effTkArr=[]
        for i in range(len(tkArr)):
            print("Thickness = %s nm; Band-gap = %s eV; Material = %s " % (str(tkArr[i]),str(bgArrReal[j]),fileName[j]))
            maxEff=calEff(tkArr[i],bgArrReal[j],bgShift[j],0,0,1,fileName[j],bgArrFun[j])
            print(maxEff)
            effTkArr.append(maxEff) 
        effDict[fileName[j]]=effTkArr 


File=pd.DataFrame(effDict,index=tkArr)
File.to_csv('MAXEFF.txt', header=True, index=True, sep=str(u'\t'), mode='w')


"""
outfile = open("test.txt", 'w')
for row in tkArr:
    outfile.write('%0.1f' % row)
    outfile.write('\n')
outfile.write('\n')
for column in effTkArr:
    outfile.write('%.13f' % column)
    outfile.write('\n')
outfile.close()
"""





"""
#n vs Eg

try:
    Material = str(sys.argv[1])
    DBgMat = float(sys.argv[2])
    IBgMat = float(sys.argv[3])
    fileName = str(sys.argv[4])
    #a = float(sys.argv[3])
    #b = float(sys.argv[4])
    #n = int(sys.argv[5])
    
except:
    print ("The Material and the bandgaps must be inputed on the command line")
    Material = str(input("Material="))
    DBgMat = float(input("Eg_dir="))
    IBgMat = float(input("Eg_ind="))
    fileName = str(input("fileName=")) 
    #a = float(input("Eg start="))
    #b = float(input("Eg end="))
    #n = int(input("points="))


bgArrReal= DBgMat
bgArrFun = bgArrReal
bgShift = IBgMat


# bandgap
Eg= np.linspace(0.4,4.0,10) 

# absorption
#fileName = "CIS_ABS.txt"

# thickness
tkArr = np.array([10,100,1000,10000])
#tkArr = np.linspace(a,b,n)


effDict={}
for j in range(len(tkArr)):
        print('------------------------------------------------------')
        print('------------------------------------------------------')
        effTkArr=[]
        for i in range(len(Eg)):
            print("Bandgap = %s eV; Thickness = %s nm; Material = %s " % (str(Eg[i]),str(tkArr[j]),fileName))
            maxEff=calEff(tkArr[j],bgArrReal,bgShift,0,0,1,fileName,bgArrFun,Eg[i])
            print(maxEff)
            effTkArr.append(maxEff)
        effDict[fileName[j]]=effTkArr

File=pd.DataFrame(effDict,index=Eg)
File.to_csv('EFFvsEg.txt', header=True, index=True, sep=str(u'\t'), mode='w')
"""


"""
# Fill factor

########### band-gap array #############
# thickness nm
#tkArr=np.array([1,2,3,4,5,6,7,8,9,10]) #small
#tkArr=np.array([5,10,25,50,100,250,500,1000,2000,5000,10000,50000,100000,500000,1000000]) #big
#tkArr = np.array([10,100,1000,10000])
tkArr = np.array([100])
        
# band gap of material eV
Eg= np.linspace(0.4,4.0,10)
#Eg = np.array([1.0])
try:
    Material = str(sys.argv[1])
    DBgMat = float(sys.argv[2])
    IBgMat = float(sys.argv[3])
    fileName = [str(sys.argv[4])]
    
except:
    print ("The Material and the bandgaps must be inputed on the command line")
    Material = str(input("Material="))
    DBgMat = float(input("Eg_dir="))
    IBgMat = float(input("Eg_ind=")) 
    fileName = [str(input("fileName="))] 

        
bgArrReal=np.array([DBgMat]*len(fileName)) # minimum direct bandgap
bgArrFun=bgArrReal
bgShift=np.array([IBgMat]*len(fileName)) # minimum bandgap (can be indirect)

effDict={}
effTkArr=[]
Voc = []
Jsc = []
FF = []

for j in range(len(bgArrReal)):
    print('------------------------------------------------------')
    print('------------------------------------------------------')
    for i in range(len(Eg)):
        print("Thickness = %s nm; Band-gap = %s eV; Material = %s " % (str(tkArr[j]),str(bgArrReal[j]),fileName[j]))
        maxEff=calEff(tkArr[j],bgArrReal[j],bgShift[j],0,0,1,fileName[j],bgArrFun[j],Eg[i])
        print(maxEff)
        Voc.append(maxEff[0])
        Jsc.append(maxEff[1])
        FF.append(maxEff[2]) 
    effTkArr = [Voc,Jsc,FF]
    effDict[fileName[j]]=effTkArr 
    



#File=pd.DataFrame(effDict,index=Voc)
#File.to_csv('FillFactor.txt', header=True, index=True, sep=str(u'\t'), mode='w')

outfile = open("test.txt", 'w')
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
"""

"""
# Complete code?

########### band-gap array #############
# thickness nm
#tkArr=np.array([1,2,3,4,5,6,7,8,9,10]) #small
#tkArr=np.array([5,10,25,50,100,250,500,1000,2000,5000,10000,50000,100000,500000,1000000]) #big
tkArr = np.array([10,100,1000,10000])
#tkArr = np.array([100])
        
# band gap of material eV
Eg= np.linspace(0.4,4.0,10)
#Eg = np.array([1.0])
try:
    Material = str(sys.argv[1])
    DBgMat = float(sys.argv[2])
    IBgMat = float(sys.argv[3])
    fileName = [str(sys.argv[4])]
    
except:
    print ("The Material and the bandgaps must be inputed on the command line")
    Material = str(input("Material="))
    DBgMat = float(input("Eg_dir="))
    IBgMat = float(input("Eg_ind=")) 
    fileName = [str(input("fileName="))] 

        
bgArrReal=np.array([DBgMat]*len(fileName)) # minimum direct bandgap
bgArrFun=bgArrReal
bgShift=np.array([IBgMat]*len(fileName)) # minimum bandgap (can be indirect)

effDict={}
effTkArr=[]
Voc = []
Jsc = []
FF = []

for j in range(len(bgArrReal)):
    print('------------------------------------------------------')
    print('------------------------------------------------------')
    for i in range(len(Eg)):
        print("Thickness = %s nm; Band-gap = %s eV; Material = %s " % (str(tkArr[j]),str(bgArrReal[j]),fileName[j]))
        maxEff=calEff(tkArr[j],bgArrReal[j],bgShift[j],0,0,1,fileName[j],bgArrFun[j],Eg[i])[0]
        print(maxEff)
        effTkArr.append(maxEff)
    effDict[fileName[j]]=effTkArr 
    



File=pd.DataFrame(effDict,index=effTkArr)
File.to_csv('FillFactor.txt', header=True, index=True, sep=str(u'\t'), mode='w')
"""

"""
outfile = open("test.txt", 'w')
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
"""
