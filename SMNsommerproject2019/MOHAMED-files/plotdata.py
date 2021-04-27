# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 13:50:44 2019

@author: mohe9
"""

from numpy import *
from matplotlib.pyplot import *
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

