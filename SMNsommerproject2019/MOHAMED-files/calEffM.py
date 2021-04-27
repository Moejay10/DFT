from __future__ import division, print_function, unicode_literals
import numpy, scipy.interpolate, scipy.integrate
from pylab import *
import matplotlib.pyplot as plt
import bigfloat
import time

import pandas as pd

import sys

"""
# Original

############# semi theoretical modelling ############
def calEff(tkCell,Bg,shiftEne,relativeRatio,relativeValue,realMaterial,fileName,gapFun):
#start_time = time.time()

# 1 eV = 1.60217662E-19 joules
    ev2j=1.60217662E-19
# 1 joule = 1/ev2j eV
    j2ev=6.24150913E18 
# elementary charge = 1.60217662E-19 (Coulombs)
    qChar=1.60217662E-19
# planck constant 6.626070040E-34 (j*s) 
    plkC=6.626070040E-34
# the speed of light 2.99792458E8 (nm/s)
    speedL=2.99792458E17
# room temperature 300 (K)
    tCell = 300
# Boltzmann constant 1.3806488E-23 (J/K)
    blmC=1.3806488E-23

# thickness of layer (nm) 
#tkCell=500

# load absorption coefficient
#data = np.loadtxt(dataFile)
# fake points (0) added in the first three lines, numerical reason
    dataAM=np.loadtxt('AM.txt') #('/Users/rongzhenchen/workplace/bookChapter/DIELECTRICFUNCTION/pd/data/AM.txt')

    # band-gap of loaded material (eV)

    #kb=0.02585...
    broaden=0.000001
    #    broaden=0.01
    #?E
    stepE= 0.001

    def sigmoid(x,a,b):
    # sigmoid function with parameters a = center; b = width
    #        return 1./(1.+bigfloat.exp(-(x-a)/b,bigfloat.precision(100)))
       return 1./(1.+exp(-(x-a)/b))
    # Wavelength is in column 0, AM1.5G data is column 2

    #    f, axarr = plt.subplots(2, sharex=True)

    AM15 = dataAM[:,[0,2]]


    shiftAM=0 #nm
    shiftAME=0 # eV
    AM15[:,0]=AM15[:,0]+shiftAM

    tempEE=1240/AM15[:,0]+shiftAME
    newWL=1240/tempEE
    AM15[:,0]=newWL


    #print(AM15[0:10,:])
    AM15interp = scipy.interpolate.interp1d(AM15[:,0], AM15[:,1])
    
    
    



    # wavelength unit nm
    wavelength_min = 280+shiftAM
    wavelength_max = 4000+shiftAM
    
    # unit in eV
    E_min = plkC*speedL/ev2j/wavelength_max+shiftAME
    E_max = plkC*speedL/ev2j/wavelength_min+shiftAME
    
    #E_min = 0.0
    #E_max = 25.0

    ############# semi theoretical modelling ############
    if(realMaterial):
        #print('real absorption')
        #print(fileName)
        data = np.loadtxt(fileName)
        #data = fileName
        lenData=len(data[:,0])
        absM=np.ndarray(shape=(lenData,2), dtype=float)
    # first column energy in eV, second column absorption coefficient in nm**(-1)
        absM[:,0] = data[:,0]
        absM[:,1] = (data[:,1])/1000. 
        absM[:,1] =  sigmoid(absM[:,0], Bg, broaden)*absM[:,1]
        
        if(shiftEne >= Bg):
            absM[:,0]=absM[:,0]+(shiftEne-Bg)
            x=np.insert(absM[:,0],0,0)
            y=np.insert(absM[:,1],0,0)
        else:
            absM[:,0]=absM[:,0]+(shiftEne-Bg)
            absM[:,0]=(absM[:,0]>=0)*absM[:,0]
            x=absM[:,0] 
            y=absM[:,1] 

        absE_max=x[x.argmax()]
    else:
    ############# pure theoretical modelling ############
        absE_max=10
        ##### (x-ctPt)**(1./2)/x, 0.001.
        x=(np.linspace(0.001,absE_max,1000001))
        absM=numpy.zeros(shape=(1000001,2))

        absM[:,0]=x

        ####### constant
        
        

        #######  (x-ctPt)**(1./2)/x, ctPt=0.01, 
        ######  CHANGE 0.000
        ctPt=0.01
        absM[:,1]=5E-4*(x>=ctPt)*(abs(x-ctPt))**(1./2)/x
        absM[:,0]=absM[:,0]+shiftEne-ctPt
        

        x=np.insert(absM[:,0],0,0)
        y=np.insert(absM[:,1],0,0)


     ######### test 
     ######### f=C*exp(-(E-Ec).^2/Ec)   
     ######### a(E)=a(E).*f




    if(relativeRatio):
        #print("relativeRatio is true")
        absE_max=relativeValue
    #else:
        #print("relativeRatio is false")

    ABSinterp = scipy.interpolate.interp1d(x, y)
    
    # A(E) = 1 or 0
    def A1(E):
      if E > Bg:
        return 1
      else:
        return 0

    # A(E) 
    def A(E):
      if E < Bg:
        return 0
      else:
        return (1 - exp(-2*ABSinterp(E)*tkCell))

    
    # A(E) 
    def A(E):
      return (1 - exp(-2*ABSinterp(E)*tkCell))

    ###################### total income light  ##############################
    def SPhotonsPerTEA(E):
        wavelength = plkC*speedL/ev2j/E # E in eV, wavelength in nm
        return AM15interp(wavelength)*(1./E)*(1./E**2)

    PowerPerTEA = lambda E : E * SPhotonsPerTEA(E)
    # solar_constant unit: W/m**2
    solar_constant =plkC*speedL*j2ev*scipy.integrate.quad(PowerPerTEA,E_min,E_max,epsabs=1e-100,epsrel=1e-100, full_output=1)[0] 


    ############# functions for simpson integral
    def Get_N(a,b,width):
    # width is step interval
        N=int((b-a)/width + 1)
        if N%2 == 0:
            N=N+1
        return N

    def GenerateData1(a,b,n,width):
        datas = []
        r=a
        for i in range(0,n):
            datas.append(func1(r))
            r = r+width
        return datas

    def GenerateData2(a,b,n,width):
        datas = []
        r=a
        for i in range(0,n):
            datas.append(func2(r))
            r = r+width
        return datas

    def simpson_integral(datas,width,n):
        sum = datas[0]+datas[n-1]
        for i in range(2,n):
            if i%2== 0:
                sum = sum +4*datas[i-1]
            else:
                sum = sum +2*datas[i-1]
        return sum*width/3.0
        
    #################### integral cut off #########
    if(relativeRatio):
    
        ######################  photons above bandgap ######
        def func1(E): 
            wavelength = plkC*speedL/ev2j/E
            return A(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            #return A1(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            
        ######################  radiative recombination ####
        def func2(E): 
    #         return ((1-exp(-2.*tkCell*ABSinterp(E)))*E**2./(exp(E/(blmC*tCell*j2ev))-1))
            return ((1-bigfloat.exp(-2.*tkCell*ABSinterp(E),bigfloat.precision(100)))*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(100))-1))
            #return A1(E)*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(100))-1)

        a=E_min # dn 
        #a=gapFun # dn 
        
        b=relativeValue

        width=stepE #step
        N=Get_N(a,b,width)
        datas1 = GenerateData1(a,b,N,width)
        prefactor0=plkC*speedL*j2ev*qChar*j2ev
        solar_constant1=(prefactor0*simpson_integral(datas1,width,N))

        datas2 = GenerateData2(a,b,N,width)
        prefactor1=qChar*2*pi/((speedL/1E9)**2*(plkC*j2ev)**3)
        RR0=(prefactor1*simpson_integral(datas2,width,N))

    #################### normal #########
    else:
        ######################  photons above bandgap ######
        def func(E): 
            wavelength = plkC*speedL/ev2j/E
            return A(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            #return A1(E)*AM15interp(wavelength)*(1./E)* (1./E**2)

        func1=lambda E : func(E) 
        prefactor0=plkC*speedL*j2ev*qChar*j2ev
        
        #solar_constant1 = prefactor0*scipy.integrate.quad(func1,gapFun,E_max,epsabs=1e-100,epsrel=1e-100,full_output=1)[0]
        solar_constant1 = prefactor0*scipy.integrate.quad(func1,E_min,E_max,epsabs=1e-100,epsrel=1e-100,full_output=1)[0] 
        #print('')
        #print('-------solar_constant1 A/m^2')
        #print(solar_constant1)

        ######################  radiative recombination ####
        def SPhotonsPerTEA3(E):
            return ((1-bigfloat.exp(-2.*tkCell*ABSinterp(E),bigfloat.precision(500)))*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(500))-1))
            #return A1(E)*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(500))-1)
    #        return prefactor1*(1-exp(-2*ABSinterp(E)*tkCell))*E**2/(exp(E/(blmC*tCell))-1)

        integrand = lambda E : SPhotonsPerTEA3(E)
        prefactor1=qChar*2*pi/((speedL/1E9)**2*(plkC*j2ev)**3)
        
        RR0 = prefactor1*scipy.integrate.quad(integrand, E_min, E_max, epsabs=1e-100,epsrel=1e-100,full_output=1)[0]
        #RR0 = prefactor1*scipy.integrate.quad(integrand, gapFun, E_max, epsabs=1e-100,epsrel=1e-100,full_output=1)[0]        
        
        na=2.51E19 # 1/cm^3
        Eg=gapFun # eV CuSbSeTe
        Cauger=1E-31 # cm^6/s
        Ndefect=1E16
        prefactor2=qChar*Cauger*na*na*tkCell*1E-7*1E4*Ndefect*exp(-Eg*ev2j/(blmC*tCell)) # A/m^2
        prefactor3=qChar*Cauger*na*na*na*tkCell*1E-7*1E4*exp(-3*Eg*ev2j/(2*blmC*tCell)) # A/m^2 
    ###################### find maximum  ##############################
    def current_density(V):
        # original
        #return (solar_constant1 - RR0*(exp(V*qChar/(blmC*tCell))))
        # paper == clas 1
        #return (solar_constant1 - RR0*(exp(V*qChar/(blmC*tCell))-1)-prefactor2*exp(3*qChar*V/(2*blmC*tCell)))
        #return (solar_constant1 - (RR0+1*prefactor2+0*prefactor3*exp(qChar*V/(2*blmC*tCell)))*exp(V*qChar/(blmC*tCell)))
        #return (solar_constant1 - ((RR0+1*prefactor2)*(exp(V*qChar/(blmC*tCell))-1)))
        return (solar_constant1 - ((RR0)*(exp(V*qChar/(blmC*tCell))-1)))
        # clas1

    from scipy.optimize import fmin 

    def fmax(func_to_maximize, initial_guess=0):
        func_to_minimize = lambda x : -func_to_maximize(x)
        return fmin(func_to_minimize, initial_guess, disp=False)[0]

    def V_mpp():
        return fmax(lambda V : V * current_density(V))
    
    
    #def J_mpp(Egap):
    # current at max power point 
    #return current_density(V_mpp())
    
  
    def max_power():
        V = V_mpp()
        return V * current_density(V)

    def max_efficiency():
        return max_power() / solar_constant

    maxEff=max_efficiency()
    return(maxEff)
    #print("--- %s seconds ---" % (time.time() - start_time))
"""



"""
n vs Eg

############# semi theoretical modelling ############
def calEff(tkCell,Bg,shiftEne,relativeRatio,relativeValue,realMaterial,fileName,gapFun,Egap):
#start_time = time.time()

# 1 eV = 1.60217662E-19 joules
    ev2j=1.60217662E-19
# 1 joule = 1/ev2j eV
    j2ev=6.24150913E18 
# elementary charge = 1.60217662E-19 (Coulombs)
    qChar=1.60217662E-19
# planck constant 6.626070040E-34 (j*s) 
    plkC=6.626070040E-34
# the speed of light 2.99792458E8 (nm/s)
    speedL=2.99792458E17
# room temperature 300 (K)
    tCell = 300
# Boltzmann constant 1.3806488E-23 (J/K)
    blmC=1.3806488E-23

# thickness of layer (nm) 
#tkCell=500

# load absorption coefficient
#data = np.loadtxt(dataFile)
# fake points (0) added in the first three lines, numerical reason
    dataAM=np.loadtxt('AM.txt') #('/Users/rongzhenchen/workplace/bookChapter/DIELECTRICFUNCTION/pd/data/AM.txt')

    # band-gap of loaded material (eV)

    #kb=0.02585...
    broaden=0.000001
    #    broaden=0.01
    #?E
    stepE=0.001

    def sigmoid(x,a,b):
    # sigmoid function with parameters a = center; b = width
    #        return 1./(1.+bigfloat.exp(-(x-a)/b,bigfloat.precision(100)))
       return 1./(1.+exp(-(x-a)/b))
    # Wavelength is in column 0, AM1.5G data is column 2

    #    f, axarr = plt.subplots(2, sharex=True)

    AM15 = dataAM[:,[0,2]]


    shiftAM=0 #nm
    shiftAME=0 # eV
    AM15[:,0]=AM15[:,0]+shiftAM

    tempEE=1240/AM15[:,0]+shiftAME
    newWL=1240/tempEE
    AM15[:,0]=newWL


    #print(AM15[0:10,:])
    AM15interp = scipy.interpolate.interp1d(AM15[:,0], AM15[:,1])
    
    
    



    # wavelength unit nm
    wavelength_min = 280+shiftAM
    wavelength_max = 4000+shiftAM
    
    # unit in eV
    E_min = plkC*speedL/ev2j/wavelength_max+shiftAME
    E_max = plkC*speedL/ev2j/wavelength_min+shiftAME
    
    #E_min = 0.0
    #E_max = 10.0

    ############# semi theoretical modelling ############
    if(realMaterial):
        #print('real absorption')
        #print(fileName)
        data = np.loadtxt(fileName)
        #data = fileName
        lenData=len(data[:,0])
        absM=np.ndarray(shape=(lenData,2), dtype=float)
    # first column energy in eV, second column absorption coefficient in nm**(-1)
        absM[:,0] = data[:,0]
        absM[:,1] = (data[:,1])/1000. 
        absM[:,1] =  sigmoid(absM[:,0], Bg, broaden)*absM[:,1]
        
        if(shiftEne >= Bg):
            absM[:,0]=absM[:,0]+(shiftEne-Bg)
            x=np.insert(absM[:,0],0,0)
            y=np.insert(absM[:,1],0,0)
        else:
            absM[:,0]=absM[:,0]+(shiftEne-Bg)
            absM[:,0]=(absM[:,0]>=0)*absM[:,0]
            x=absM[:,0] 
            y=absM[:,1] 

        absE_max=x[x.argmax()]
    else:
    ############# pure theoretical modelling ############
        absE_max=10
        ##### (x-ctPt)**(1./2)/x, 0.001.
        x=(np.linspace(0.001,absE_max,1000001))
        absM=numpy.zeros(shape=(1000001,2))

        absM[:,0]=x

        ####### constant
        
        

        #######  (x-ctPt)**(1./2)/x, ctPt=0.01, 
        ######  CHANGE 0.000
        ctPt=0.01
        absM[:,1]=5E-4*(x>=ctPt)*(abs(x-ctPt))**(1./2)/x
        absM[:,0]=absM[:,0]+shiftEne-ctPt
        

        x=np.insert(absM[:,0],0,0)
        y=np.insert(absM[:,1],0,0)


     ######### test 
     ######### f=C*exp(-(E-Ec).^2/Ec)   
     ######### a(E)=a(E).*f




    if(relativeRatio):
        #print("relativeRatio is true")
        absE_max=relativeValue
    #else:
        #print("relativeRatio is false")

    ABSinterp = scipy.interpolate.interp1d(x, y)
    
    # A(E) = 1 or 0, to test what makes the efficiency go down
    def A1(E):
      if E > Bg:
        return 1
      else:
        return 0
      
    # A(E) 
    def A(E):
      return (1 - exp(-2*ABSinterp(E)*tkCell))
      
    
    def Asm(E):
      return 2*tkCell*ABSinterp(E) 
    

    ###################### total income light  ##############################
    def SPhotonsPerTEA(E):
        wavelength = plkC*speedL/ev2j/E # E in eV, wavelength in nm
        return AM15interp(wavelength)*(1./E)*(1./E**2)

    PowerPerTEA = lambda E : E * SPhotonsPerTEA(E)
    # solar_constant unit: W/m**2
    solar_constant =plkC*speedL*j2ev*scipy.integrate.quad(PowerPerTEA,E_min,E_max,epsabs=1e-100,epsrel=1e-100, full_output=1)[0] 


    ############# functions for simpson integral
    def Get_N(a,b,width):
    # width is step interval
        N=int((b-a)/width + 1)
        if N%2 == 0:
            N=N+1
        return N

    def GenerateData1(a,b,n,width):
        datas = []
        r=a
        for i in range(0,n):
            datas.append(func1(r))
            r = r+width
        return datas

    def GenerateData2(a,b,n,width):
        datas = []
        r=a
        for i in range(0,n):
            datas.append(func2(r))
            r = r+width
        return datas

    def simpson_integral(datas,width,n):
        sum = datas[0]+datas[n-1]
        for i in range(2,n):
            if i%2== 0:
                sum = sum +4*datas[i-1]
            else:
                sum = sum +2*datas[i-1]
        return sum*width/3.0
        
    #################### integral cut off #########
    if(relativeRatio):
    
        ######################  photons above bandgap ######
        def func1(E): 
            wavelength = plkC*speedL/ev2j/E
            return A(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            #return Asm(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            #return A1(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            
        ######################  radiative recombination ####
        def func2(E): 
    #         return ((1-exp(-2.*tkCell*ABSinterp(E)))*E**2./(exp(E/(blmC*tCell*j2ev))-1))
            return ((1-bigfloat.exp(-2.*tkCell*ABSinterp(E),bigfloat.precision(100)))*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(100))-1))
            #return Asm(E)*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(100))-1)

        a=E_min # dn 
        #a=gapFun # dn 
        
        b=relativeValue

        width=stepE #step
        N=Get_N(a,b,width)
        datas1 = GenerateData1(a,b,N,width)
        prefactor0=plkC*speedL*j2ev*qChar*j2ev
        solar_constant1=(prefactor0*simpson_integral(datas1,width,N))

        datas2 = GenerateData2(a,b,N,width)
        prefactor1=qChar*2*pi/((speedL/1E9)**2*(plkC*j2ev)**3)
        RR0=(prefactor1*simpson_integral(datas2,width,N))

    #################### normal #########
    else:
        ######################  photons above bandgap ######
        def func(E): 
            wavelength = plkC*speedL/ev2j/E
            return A(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            #return Asm(E)*AM15interp(wavelength)*(1./E)* (1./E**2)

        func1=lambda E : func(E) 
        prefactor0=plkC*speedL*j2ev*qChar*j2ev
        
        def solar_constant1(Egap):
            return prefactor0*scipy.integrate.quad(func1,Egap,E_max,epsabs=1e-100,epsrel=1e-100,full_output=1)[0]
        
        #solar_constant1 = prefactor0*scipy.integrate.quad(func1,E_min,E_max,epsabs=1e-100,epsrel=1e-100,full_output=1)[0] 
        #print('')
        #print('-------solar_constant1 A/m^2')
        #print(solar_constant1)

        ######################  radiative recombination ####
        def SPhotonsPerTEA3(E):
            return ((1-bigfloat.exp(-2.*tkCell*ABSinterp(E),bigfloat.precision(500)))*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(500))-1))
            #return Asm(E)*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(500))-1)
    #        return prefactor1*(1-exp(-2*ABSinterp(E)*tkCell))*E**2/(exp(E/(blmC*tCell))-1)

        integrand = lambda E : SPhotonsPerTEA3(E)
        prefactor1=qChar*2*pi/((speedL/1E9)**2*(plkC*j2ev)**3)
        
        def RR0(Egap):
            return prefactor1*scipy.integrate.quad(integrand, Egap, E_max, epsabs=1e-100,epsrel=1e-100,full_output=1)[0]
        #RR0 = prefactor1*scipy.integrate.quad(integrand, E_min, E_max, epsabs=1e-100,epsrel=1e-100,full_output=1)[0]        
        
        na=2.51E19 # 1/cm^3
        Eg=gapFun # eV CuSbSeTe
        Cauger=1E-31 # cm^6/s
        Ndefect=1E16
        prefactor2=qChar*Cauger*na*na*tkCell*1E-7*1E4*Ndefect*exp(-Eg*ev2j/(blmC*tCell)) # A/m^2
        prefactor3=qChar*Cauger*na*na*na*tkCell*1E-7*1E4*exp(-3*Eg*ev2j/(2*blmC*tCell)) # A/m^2 
    ###################### find maximum  ##############################
    def current_density(V,Egap):
        # original
        #return (solar_constant1 - RR0*(exp(V*qChar/(blmC*tCell))))
        # paper == clas 1
        #return (solar_constant1 - RR0*(exp(V*qChar/(blmC*tCell))-1)-prefactor2*exp(3*qChar*V/(2*blmC*tCell)))
        #return (solar_constant1 - (RR0+1*prefactor2+0*prefactor3*exp(qChar*V/(2*blmC*tCell)))*exp(V*qChar/(blmC*tCell)))
        #return (solar_constant1 - ((RR0+1*prefactor2)*(exp(V*qChar/(blmC*tCell))-1)))
        return (solar_constant1(Egap) - ((RR0(Egap))*(exp(V*qChar/(blmC*tCell))-1)))
        # clas1

    from scipy.optimize import fmin 

    def fmax(func_to_maximize, initial_guess=0):
        func_to_minimize = lambda x : -func_to_maximize(x)
        return fmin(func_to_minimize, initial_guess, disp=False)[0]

    def V_mpp(Egap):
        return fmax(lambda V : V * current_density(V,Egap))
    
    
    def J_mpp(Egap):
    # current at max power point 
        return current_density(V_mpp(Egap),Egap)
    
    
    
  
    def max_power(Egap):
        V = V_mpp(Egap)
        return V * current_density(V,Egap)

    def max_efficiency(Egap):
        return max_power(Egap) / solar_constant

    maxEff=max_efficiency(Egap)
    return(maxEff)
    #print("--- %s seconds ---" % (time.time() - start_time))
"""


# Fill factor

############# semi theoretical modelling ############
def calEff(tkCell,Bg,shiftEne,relativeRatio,relativeValue,realMaterial,fileName,gapFun,Egap):
#start_time = time.time()

# 1 eV = 1.60217662E-19 joules
    ev2j=1.60217662E-19
# 1 joule = 1/ev2j eV
    j2ev=6.24150913E18 
# elementary charge = 1.60217662E-19 (Coulombs)
    qChar=1.60217662E-19
# planck constant 6.626070040E-34 (j*s) 
    plkC=6.626070040E-34
# the speed of light 2.99792458E8 (nm/s)
    speedL=2.99792458E17
# room temperature 300 (K)
    tCell = 300
# Boltzmann constant 1.3806488E-23 (J/K)
    blmC=1.3806488E-23

# thickness of layer (nm) 
#tkCell=500

# load absorption coefficient
#data = np.loadtxt(dataFile)
# fake points (0) added in the first three lines, numerical reason
    dataAM=np.loadtxt('AM.txt') #('/Users/rongzhenchen/workplace/bookChapter/DIELECTRICFUNCTION/pd/data/AM.txt')

    # band-gap of loaded material (eV)

    #kb=0.02585...
    broaden=0.000001
    #    broaden=0.01
    #?E
    stepE=0.001

    def sigmoid(x,a,b):
    # sigmoid function with parameters a = center; b = width
    #        return 1./(1.+bigfloat.exp(-(x-a)/b,bigfloat.precision(100)))
       return 1./(1.+exp(-(x-a)/b))
    # Wavelength is in column 0, AM1.5G data is column 2

    #    f, axarr = plt.subplots(2, sharex=True)

    AM15 = dataAM[:,[0,2]]


    shiftAM=0 #nm
    shiftAME=0 # eV
    AM15[:,0]=AM15[:,0]+shiftAM

    tempEE=1240/AM15[:,0]+shiftAME
    newWL=1240/tempEE
    AM15[:,0]=newWL


    #print(AM15[0:10,:])
    AM15interp = scipy.interpolate.interp1d(AM15[:,0], AM15[:,1])
    
    
    



    # wavelength unit nm
    wavelength_min = 280+shiftAM
    wavelength_max = 4000+shiftAM
    
    # unit in eV
    E_min = plkC*speedL/ev2j/wavelength_max+shiftAME
    #E_max = plkC*speedL/ev2j/wavelength_min+shiftAME
    
    #E_min = 0.0
    E_max = 10.0

    ############# semi theoretical modelling ############
    if(realMaterial):
        #print('real absorption')
        #print(fileName)
        data = np.loadtxt(fileName)
        #data = fileName
        lenData=len(data[:,0])
        absM=np.ndarray(shape=(lenData,2), dtype=float)
    # first column energy in eV, second column absorption coefficient in nm**(-1)
        absM[:,0] = data[:,0]
        absM[:,1] = (data[:,1])/1000. 
        #absM[:,1] =  sigmoid(absM[:,0], Bg, broaden)*absM[:,1]
        
        if(shiftEne >= Bg):
            absM[:,0]=absM[:,0]+(shiftEne-Bg)
            x=np.insert(absM[:,0],0,0)
            y=np.insert(absM[:,1],0,0)
        else:
            absM[:,0]=absM[:,0]+(shiftEne-Bg)
            absM[:,0]=(absM[:,0]>=0)*absM[:,0]
            x=absM[:,0] 
            y=absM[:,1] 

        absE_max=x[x.argmax()]
    else:
    ############# pure theoretical modelling ############
        absE_max=10
        ##### (x-ctPt)**(1./2)/x, 0.001.
        x=(np.linspace(0.001,absE_max,1000001))
        absM=numpy.zeros(shape=(1000001,2))

        absM[:,0]=x

        ####### constant
        
        

        #######  (x-ctPt)**(1./2)/x, ctPt=0.01, 
        ######  CHANGE 0.000
        ctPt=0.01
        absM[:,1]=5E-4*(x>=ctPt)*(abs(x-ctPt))**(1./2)/x
        absM[:,0]=absM[:,0]+shiftEne-ctPt
        

        x=np.insert(absM[:,0],0,0)
        y=np.insert(absM[:,1],0,0)


     ######### test 
     ######### f=C*exp(-(E-Ec).^2/Ec)   
     ######### a(E)=a(E).*f




    if(relativeRatio):
        #print("relativeRatio is true")
        absE_max=relativeValue
    #else:
        #print("relativeRatio is false")

    ABSinterp = scipy.interpolate.interp1d(x, y)
    
    # A(E) = 1 or 0, to test what makes the efficiency go down
    def A1(E):
      if E > Bg:
        return 1
      else:
        return 0
      
     # A(E) 
    def A(E):
      if E < Bg:
        return 0
      else:
        #return (1 - exp(-2*ABSinterp(E)*tkCell))
        return (1-bigfloat.exp(-2.*tkCell*ABSinterp(E),bigfloat.precision(500)))
      
    
    def Asm(E):
      return 2*tkCell*ABSinterp(E) 
    

    ###################### total income light  ##############################
    def SPhotonsPerTEA(E):
        wavelength = plkC*speedL/ev2j/E # E in eV, wavelength in nm
        return AM15interp(wavelength)*(1./E)*(1./E**2)

    PowerPerTEA = lambda E : E * SPhotonsPerTEA(E)
    # solar_constant unit: W/m**2
    solar_constant =plkC*speedL*j2ev*scipy.integrate.quad(PowerPerTEA,E_min,E_max,epsabs=1e-100,epsrel=1e-100, full_output=1)[0] 


    ############# functions for simpson integral
    def Get_N(a,b,width):
    # width is step interval
        N=int((b-a)/width + 1)
        if N%2 == 0:
            N=N+1
        return N

    def GenerateData1(a,b,n,width):
        datas = []
        r=a
        for i in range(0,n):
            datas.append(func1(r))
            r = r+width
        return datas

    def GenerateData2(a,b,n,width):
        datas = []
        r=a
        for i in range(0,n):
            datas.append(func2(r))
            r = r+width
        return datas

    def simpson_integral(datas,width,n):
        sum = datas[0]+datas[n-1]
        for i in range(2,n):
            if i%2== 0:
                sum = sum +4*datas[i-1]
            else:
                sum = sum +2*datas[i-1]
        return sum*width/3.0
        
    #################### integral cut off #########
    if(relativeRatio):
    
        ######################  photons above bandgap ######
        def func1(E): 
            wavelength = plkC*speedL/ev2j/E
            return A(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            #return Asm(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            #return A1(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            
        ######################  radiative recombination ####
        def func2(E): 
    #         return ((1-exp(-2.*tkCell*ABSinterp(E)))*E**2./(exp(E/(blmC*tCell*j2ev))-1))
            return ((1-bigfloat.exp(-2.*tkCell*ABSinterp(E),bigfloat.precision(100)))*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(100))-1))
            #return Asm(E)*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(100))-1)

        a=E_min # dn 
        #a=gapFun # dn 
        
        b=relativeValue

        width=stepE #step
        N=Get_N(a,b,width)
        datas1 = GenerateData1(a,b,N,width)
        prefactor0=plkC*speedL*j2ev*qChar*j2ev
        solar_constant1=(prefactor0*simpson_integral(datas1,width,N))

        datas2 = GenerateData2(a,b,N,width)
        prefactor1=qChar*2*pi/((speedL/1E9)**2*(plkC*j2ev)**3)
        RR0=(prefactor1*simpson_integral(datas2,width,N))

    #################### normal #########
    else:
        ######################  photons above bandgap ######
        def func(E): 
            wavelength = plkC*speedL/ev2j/E
            return A(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            #return Asm(E)*AM15interp(wavelength)*(1./E)* (1./E**2)

        func1=lambda E : func(E) 
        prefactor0=plkC*speedL*j2ev*qChar*j2ev
        
        def solar_constant1(Egap):
            return prefactor0*scipy.integrate.quad(func1,Egap,E_max,epsabs=1e-100,epsrel=1e-10,full_output=1)[0]
        
        #solar_constant1 = prefactor0*scipy.integrate.quad(func1,E_min,E_max,epsabs=1e-100,epsrel=1e-100,full_output=1)[0] 
        #print('')
        #print('-------solar_constant1 A/m^2')
        #print(solar_constant1)

        ######################  radiative recombination ####
        def SPhotonsPerTEA3(E):
            return ((1-bigfloat.exp(-2.*tkCell*ABSinterp(E),bigfloat.precision(500)))*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(500))-1))
            #return Asm(E)*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(500))-1)
    #        return prefactor1*(1-exp(-2*ABSinterp(E)*tkCell))*E**2/(exp(E/(blmC*tCell))-1)

        integrand = lambda E : SPhotonsPerTEA3(E)
        prefactor1=qChar*2*pi/((speedL/1E9)**2*(plkC*j2ev)**3)
        
        def RR0(Egap):
            return prefactor1*scipy.integrate.quad(integrand, Egap, E_max, epsabs=1e-100,epsrel=1e-10,full_output=1)[0]
        #RR0 = prefactor1*scipy.integrate.quad(integrand, E_min, E_max, epsabs=1e-100,epsrel=1e-100,full_output=1)[0]        
        
        na=2.51E19 # 1/cm^3
        Eg=gapFun # eV CuSbSeTe
        Cauger=1E-31 # cm^6/s
        Ndefect=1E16
        prefactor2=qChar*Cauger*na*na*tkCell*1E-7*1E4*Ndefect*exp(-Eg*ev2j/(blmC*tCell)) # A/m^2
        prefactor3=qChar*Cauger*na*na*na*tkCell*1E-7*1E4*exp(-3*Eg*ev2j/(2*blmC*tCell)) # A/m^2 
    ###################### find maximum  ##############################
    def current_density(V,Egap):
        # original
        #return (solar_constant1 - RR0*(exp(V*qChar/(blmC*tCell))))
        # paper == clas 1
        #return (solar_constant1 - RR0*(exp(V*qChar/(blmC*tCell))-1)-prefactor2*exp(3*qChar*V/(2*blmC*tCell)))
        #return (solar_constant1 - (RR0+1*prefactor2+0*prefactor3*exp(qChar*V/(2*blmC*tCell)))*exp(V*qChar/(blmC*tCell)))
        #return (solar_constant1 - ((RR0+1*prefactor2)*(exp(V*qChar/(blmC*tCell))-1)))
        return (solar_constant1(Egap) - ((RR0(Egap))*(exp(V*qChar/(blmC*tCell))-1)))
        # clas1

    from scipy.optimize import fmin 

    def fmax(func_to_maximize, initial_guess=0):
        func_to_minimize = lambda x : -func_to_maximize(x)
        return fmin(func_to_minimize, initial_guess, disp=False)[0]

    def V_mpp(Egap):
        return fmax(lambda V : V * current_density(V,Egap))
    
    
    def J_mpp(Egap):
    # current at max power point 
        return current_density(V_mpp(Egap),Egap)
    
    
    def JSC(Egap):
        return current_density(0,Egap)
    
    def VOC(Egap):
        return (blmC*tCell/qChar)*log(solar_constant1(Egap)/RR0(Egap))
    
  
    def max_power(Egap):
        V = V_mpp(Egap)
        return V * current_density(V,Egap)
    
    def fill_factor(Egap):
        return max_power(Egap)/(JSC(Egap)*VOC(Egap))

    def max_efficiency(Egap):
        return max_power(Egap) / solar_constant

    maxEff=max_efficiency(Egap)
    
    Jsc = JSC(Egap)
    Voc = VOC(Egap)
    FF = fill_factor(Egap)
    
    return (maxEff,Voc,Jsc,FF)
    #print("--- %s seconds ---" % (time.time() - start_time))



"""
# Fixing n

############# semi theoretical modelling ############
def calEff(tkCell,Bg,shiftEne,relativeRatio,relativeValue,realMaterial,fileName,gapFun):
#start_time = time.time()

# 1 eV = 1.60217662E-19 joules
    ev2j=1.60217662E-19
# 1 joule = 1/ev2j eV
    j2ev=6.24150913E18 
# elementary charge = 1.60217662E-19 (Coulombs)
    qChar=1.60217662E-19
# planck constant 6.626070040E-34 (j*s) 
    plkC=6.626070040E-34
# the speed of light 2.99792458E8 (nm/s)
    speedL=2.99792458E17
# room temperature 300 (K)
    tCell = 300
# Boltzmann constant 1.3806488E-23 (J/K)
    blmC=1.3806488E-23

# thickness of layer (nm) 
#tkCell=500

# load absorption coefficient
#data = np.loadtxt(dataFile)
# fake points (0) added in the first three lines, numerical reason
    dataAM=np.loadtxt('AM.txt') #('/Users/rongzhenchen/workplace/bookChapter/DIELECTRICFUNCTION/pd/data/AM.txt')

    # band-gap of loaded material (eV)

    #kb=0.02585...
    broaden=0.000001
    #    broaden=0.01
    #?E
    stepE= 0.001

    def sigmoid(x,a,b):
    # sigmoid function with parameters a = center; b = width
    #        return 1./(1.+bigfloat.exp(-(x-a)/b,bigfloat.precision(100)))
       return 1./(1.+exp(-(x-a)/b))
    # Wavelength is in column 0, AM1.5G data is column 2

    #    f, axarr = plt.subplots(2, sharex=True)

    AM15 = dataAM[:,[0,2]]


    shiftAM=0 #nm
    shiftAME=0 # eV
    AM15[:,0]=AM15[:,0]+shiftAM

    tempEE=1240/AM15[:,0]+shiftAME
    newWL=1240/tempEE
    AM15[:,0]=newWL


    #print(AM15[0:10,:])
    AM15interp = scipy.interpolate.interp1d(AM15[:,0], AM15[:,1])
    
    
    



    # wavelength unit nm
    wavelength_min = 280+shiftAM
    wavelength_max = 4000+shiftAM
    
    # unit in eV
    E_min = plkC*speedL/ev2j/wavelength_max+shiftAME
    #E_max = plkC*speedL/ev2j/wavelength_min+shiftAME
    
    #E_min = 0.0
    E_max = 10.0

    ############# semi theoretical modelling ############
    if(realMaterial):
        #print('real absorption')
        #print(fileName)
        data = np.loadtxt(fileName)
        #data = fileName
        lenData=len(data[:,0])
        absM=np.ndarray(shape=(lenData,2), dtype=float)
    # first column energy in eV, second column absorption coefficient in nm**(-1)
        absM[:,0] = data[:,0]
        absM[:,1] = (data[:,1])/1000. 
        #absM[:,1] =  sigmoid(absM[:,0], Bg, broaden)*absM[:,1]
        
        if(shiftEne >= Bg):
            absM[:,0]=absM[:,0]+(shiftEne-Bg)
            x=np.insert(absM[:,0],0,0)
            y=np.insert(absM[:,1],0,0)
        else:
            absM[:,0]=absM[:,0]+(shiftEne-Bg)
            absM[:,0]=(absM[:,0]>=0)*absM[:,0]
            x=absM[:,0] 
            y=absM[:,1] 

        absE_max=x[x.argmax()]
    else:
    ############# pure theoretical modelling ############
        absE_max=10
        ##### (x-ctPt)**(1./2)/x, 0.001.
        x=(np.linspace(0.001,absE_max,1000001))
        absM=numpy.zeros(shape=(1000001,2))

        absM[:,0]=x

        ####### constant
        
        

        #######  (x-ctPt)**(1./2)/x, ctPt=0.01, 
        ######  CHANGE 0.000
        ctPt=0.01
        absM[:,1]=5E-4*(x>=ctPt)*(abs(x-ctPt))**(1./2)/x
        absM[:,0]=absM[:,0]+shiftEne-ctPt
        

        x=np.insert(absM[:,0],0,0)
        y=np.insert(absM[:,1],0,0)


     ######### test 
     ######### f=C*exp(-(E-Ec).^2/Ec)   
     ######### a(E)=a(E).*f




    if(relativeRatio):
        #print("relativeRatio is true")
        absE_max=relativeValue
    #else:
        #print("relativeRatio is false")

    ABSinterp = scipy.interpolate.interp1d(x, y)
    
    # A(E) = 1 or 0
    def A1(E):
      if E > Bg:
        return 1
      else:
        return 0
    
    # A(E) 
    def A(E):
      if E < Bg:
        return 0
      else:
        #return (1 - exp(-2*ABSinterp(E)*tkCell))
        return (1-bigfloat.exp(-2.*tkCell*ABSinterp(E),bigfloat.precision(500)))
    
    
    
    # A(E) 
    def A(E):
      return (1 - exp(-2*ABSinterp(E)*tkCell))
    
    ###################### total income light  ##############################
    def SPhotonsPerTEA(E):
        wavelength = plkC*speedL*j2ev/E # E in eV, wavelength in nm
        return AM15interp(wavelength)*(1./E)*(1.*plkC*speedL*j2ev/E**2)

    PowerPerTEA = lambda E : E * SPhotonsPerTEA(E)
    # solar_constant unit: W/m**2
    solar_constant = scipy.integrate.quad(PowerPerTEA,E_min,E_max,epsabs=1e-100,epsrel=1e-10, full_output=1)[0] 


    ############# functions for simpson integral
    def Get_N(a,b,width):
    # width is step interval
        N=int((b-a)/width + 1)
        if N%2 == 0:
            N=N+1
        return N

    def GenerateData1(a,b,n,width):
        datas = []
        r=a
        for i in range(0,n):
            datas.append(func1(r))
            r = r+width
        return datas

    def GenerateData2(a,b,n,width):
        datas = []
        r=a
        for i in range(0,n):
            datas.append(func2(r))
            r = r+width
        return datas

    def simpson_integral(datas,width,n):
        sum = datas[0]+datas[n-1]
        for i in range(2,n):
            if i%2== 0:
                sum = sum +4*datas[i-1]
            else:
                sum = sum +2*datas[i-1]
        return sum*width/3.0
        
    #################### integral cut off #########
    if(relativeRatio):
    
        ######################  photons above bandgap ######
        def func1(E): 
            wavelength = plkC*speedL/ev2j/E
            return A(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            #return A1(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            
        ######################  radiative recombination ####
        def func2(E): 
    #         return ((1-exp(-2.*tkCell*ABSinterp(E)))*E**2./(exp(E/(blmC*tCell*j2ev))-1))
            return ((1-bigfloat.exp(-2.*tkCell*ABSinterp(E),bigfloat.precision(100)))*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(100))-1))
            #return A1(E)*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(100))-1)

        a=E_min # dn 
        #a=gapFun # dn 
        
        b=relativeValue

        width=stepE #step
        N=Get_N(a,b,width)
        datas1 = GenerateData1(a,b,N,width)
        prefactor0=plkC*speedL*j2ev*qChar*j2ev
        solar_constant1=(prefactor0*simpson_integral(datas1,width,N))

        datas2 = GenerateData2(a,b,N,width)
        prefactor1=qChar*2*pi/((speedL/1E9)**2*(plkC*j2ev)**3)
        RR0=(prefactor1*simpson_integral(datas2,width,N))

    #################### normal #########
    else:
        ######################  photons above bandgap ######
        def func(E): 
            wavelength = plkC*speedL*j2ev/E
            return A(E)*AM15interp(wavelength)*(1./E)* (1./E**2)
            #return A1(E)*AM15interp(wavelength)*(1./E)* (1./E**2)

        func1=lambda E : func(E) 
        prefactor0=plkC*speedL*j2ev*qChar*j2ev
        #prefactor0=qChar*2*pi/((speedL/1E9)**2*(plkC*j2ev)**3)
        
        #solar_constant1 = prefactor0*scipy.integrate.quad(func1,gapFun,E_max,epsabs=1e-100,epsrel=1e-100,full_output=1)[0]
        solar_constant1 = prefactor0*scipy.integrate.quad(func1,E_min,E_max,epsabs=1e-100,epsrel=1e-10,full_output=1)[0] 
        #print('')
        #print('-------solar_constant1 A/m^2')
        #print(solar_constant1)

        ######################  radiative recombination ####
        def SPhotonsPerTEA3(E):
            #return ((1-bigfloat.exp(-2.*tkCell*ABSinterp(E),bigfloat.precision(500)))*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(500))-1))
            #return A(E)*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(500))-1)
            return A(E)*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(500))-1)
            #return A(E)*E**2./(exp(E/(blmC*tCell))-1)
    #        return prefactor1*(1-exp(-2*ABSinterp(E)*tkCell))*E**2/(exp(E/(blmC*tCell))-1)

        integrand = lambda E : SPhotonsPerTEA3(E)
        prefactor1=qChar*2*pi/((speedL/1E9)**2*(plkC*j2ev)**3)
        #prefactor1=qChar*2*pi/((speedL)**2*(plkC*j2ev)**3)
        
        RR0 = prefactor1*scipy.integrate.quad(integrand, E_min, E_max, epsabs=1e-100,epsrel=1e-10,full_output=1)[0]
        #RR0 = prefactor1*scipy.integrate.quad(integrand, gapFun, E_max, epsabs=1e-100,epsrel=1e-100,full_output=1)[0]        
        
        na=2.51E19 # 1/cm^3
        Eg=gapFun # eV CuSbSeTe
        Cauger=1E-31 # cm^6/s
        Ndefect=1E16
        prefactor2=qChar*Cauger*na*na*tkCell*1E-7*1E4*Ndefect*exp(-Eg*ev2j/(blmC*tCell)) # A/m^2
        prefactor3=qChar*Cauger*na*na*na*tkCell*1E-7*1E4*exp(-3*Eg*ev2j/(2*blmC*tCell)) # A/m^2 
    ###################### find maximum  ##############################
    def current_density(V):
        # original
        #return (solar_constant1 - RR0*(exp(V*qChar/(blmC*tCell))))
        # paper == clas 1
        #return (solar_constant1 - RR0*(exp(V*qChar/(blmC*tCell))-1)-prefactor2*exp(3*qChar*V/(2*blmC*tCell)))
        #return (solar_constant1 - (RR0+1*prefactor2+0*prefactor3*exp(qChar*V/(2*blmC*tCell)))*exp(V*qChar/(blmC*tCell)))
        #return (solar_constant1 - ((RR0+1*prefactor2)*(exp(V*qChar/(blmC*tCell))-1)))
        return (solar_constant1 - ((RR0)*(exp(V*qChar/(blmC*tCell))-1)))
        # clas1

    from scipy.optimize import fmin 

    def fmax(func_to_maximize, initial_guess=0):
        func_to_minimize = lambda x : -func_to_maximize(x)
        return fmin(func_to_minimize, initial_guess, disp=False)[0]

    def V_mpp():
        return fmax(lambda V : V * current_density(V))
    
    
    #def J_mpp(Egap):
    # current at max power point 
    #return current_density(V_mpp())
    
  
    def max_power():
        V = V_mpp()
        return V * current_density(V)

    def max_efficiency():
        return max_power() / solar_constant

    maxEff=max_efficiency()
    return(maxEff)
    #print("--- %s seconds ---" % (time.time() - start_time))
"""
