# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 13:50:44 2019

@author: mohe9
"""

import numpy, scipy.interpolate, scipy.integrate
import matplotlib.pyplot as plt
import sys

"""
def read(filename):
    infile = open(filename,"r")
    t_values = []
    line = infile.readline()
    words = line.split()
    v0 = float(words[-1])
    infile.readline()
    for line in infile:
        words = line.split()
        for word in words:
            t_values.append(float(word))
    infile.close()
    return v0, t_values
"""

"""
def read(filename):
    TK_values = []
    Eff_values = []
    infile = open(filename,"r")
    line = infile.readline()
    words = line.split()
    TK_values.append(float(words[-1]))
    for line in infile:
        words = line.split()
        for word in words:
            if float(word) >= 1.0:
                TK_values.append(float(word))
            else:
                Eff_values.append(float(word)) 
    infile.close()
    #return TK_values, Eff_values
    plot(TK_values,Eff_values)
    show()

read("test.txt")
"""


def read(filename):
    TK_values = []
    Eff_values = []
    infile = open(filename,"r")
    line = infile.readline()
    for line in infile:
        words = line.split()
        TK_values.append(float(words[0]))
        Eff_values.append(float(words[1])) 
    infile.close()
    plot(TK_values,Eff_values,"ro")
    show()

read("MAXEFF.txt")


"""
def sigmoid(x,a,b):
    # sigmoid function with parameters a = center; b = width
    #        return 1./(1.+bigfloat.exp(-(x-a)/b,bigfloat.precision(100)))
    return 1./(1.+exp(-(x-a)/b))

def a(fileName,Bg,shiftEne):
    
    
    #print('real absorption')
    #print(fileName)
    data = np.loadtxt(fileName)
    #data = fileName
    lenData=len(data[:,0])
    absM=np.ndarray(shape=(lenData,2), dtype=float)
    # first column energy in eV, second column absorption coefficient in nm**(-1)
    absM[:,0] = data[:,0]
    absM[:,1] = (data[:,1])/1000. 
    absM[:,1] =  sigmoid(absM[:,0], Bg, 0.001)*absM[:,1]
    
    if(shiftEne >= Bg):
            absM[:,0]=absM[:,0]+(shiftEne-Bg)
            x=np.insert(absM[:,0],0,0)
            y=np.insert(absM[:,1],0,0)
    else:
            absM[:,0]=absM[:,0]+(shiftEne-Bg)
            absM[:,0]=(absM[:,0]>=0)*absM[:,0]
            x=absM[:,0] 
            y=absM[:,1] 
    
    ABSinterp = scipy.interpolate.interp1d(absM[:,0], absM[:,1])
    
    
    E = np.linspace((1240.0/4000.0),25,1000)
    d = 100
    
    def A(E):
         return (1- exp(-2*d*ABSinterp(E)))
    
    f = scipy.integrate.quad(A,0,E[-1],epsabs=1e-100,epsrel=1e-100, full_output=1)[0]
    
    print (f)
    #plt.plot(E,ABSinterp(E))
    #plt.plot(E,f)
    #plt.show()
    
print (a("CISe_ABS.txt", 0.8495,0.8495))
print(a("CIS_ABS.txt", 0.87,0.87))
print(a("GaAs_ABS.txt",1.441,1.441))
print(a("CBS_ABS.txt", 1.61,1.35))
print(a("CGS_ABS.txt",2.0,2.0))
"""

"""
def readfile(fileName):
    
    data = np.loadtxt(fileName)
    lenData=len(data[:,0])
    absM=np.ndarray(shape=(lenData,5), dtype=float)
    # first column energy in eV, second-fifth column is efficiency 
    absM[:,0] = data[:,0]
    absM[:,1] = data[:,1]
    absM[:,1] = data[:,2]
    absM[:,1] = data[:,3]
    absM[:,1] = data[:,4]
    
    plt.plot(absM[:,0],absM[:,1])
    plt.show()

readfile("EFFvsEg.txt")
""" 

"""
# EffvsEg
def read(filename):
    Eg = []
    Eff_tk10 = []
    Eff_tk100 = []
    Eff_tk1000 = []
    Eff_tk10000 = []
    infile = open(filename,"r")
    line = infile.readline()
    #infile.readline()
    for line in infile:
        words = line.split()
        Eg.append(float(words[0]))
        Eff_tk10.append(float(words[1]))
        Eff_tk100.append(float(words[2]))
        Eff_tk1000.append(float(words[3]))
        Eff_tk10000.append(float(words[4]))
    infile.close()
    plot(Eg,Eff_tk10,"ro")
    plot(Eg,Eff_tk100,"bo")
    plot(Eg,Eff_tk1000,"go")
    plot(Eg,Eff_tk10000,"yo")
    show()

#read("EFFvsEg.txt")
"""

    