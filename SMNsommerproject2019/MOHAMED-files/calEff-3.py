from __future__ import division, print_function, unicode_literals
import numpy, scipy.interpolate, scipy.integrate
from pylab import *
import matplotlib.pyplot as plt
import bigfloat
import time

import pandas as pd 

import sys

############# semi theoretical modelling ############
def calEff(tkCell,Bg,shiftEne,relativeRatio,relativeValue,realMaterial,fileName,gapFun):


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
    dataAM=np.loadtxt('/Users/rongzhenchen/workplace/bookChapter/DIELECTRICFUNCTION/pd/data/AM.txt')

    # band-gap of loaded material (eV)

    #kb=0.02585...
    broaden=0.000001
    #    broaden=0.01
    stepE=0.001

# This function is to make the absorption coefficient data nicer. Smedge the edges of the graph.
    def sigmoid(x,a,b):
    # sigmoid function with parameters a = center; b = width
    #        return 1./(1.+bigfloat.exp(-(x-a)/b,bigfloat.precision(100)))
       return 1./(1.+exp(-(x-a)/b))
    
# Wavelength is in column 0, AM1.5G data is column 2
    AM15 = dataAM[:,[0,2]]


    shiftAM=0 #nm, it is to assign them the units I think?
    shiftAME=0 # eV, it is to assign them the units I think?
    AM15[:,0]=AM15[:,0]+shiftAM

    tempEE=1240/AM15[:,0]+shiftAME
    newWL=1240/tempEE
    AM15[:,0]=newWL # They are the same, even before the changes.


    #print(AM15[0:10,:])
    AM15interp = scipy.interpolate.interp1d(AM15[:,0], AM15[:,1])



    # wavelength unit nm
    wavelength_min = 280+shiftAM
    wavelength_max = 4000+shiftAM
    # unit in eV
    E_min = plkC*speedL/ev2j/wavelength_max+shiftAME
    E_max = plkC*speedL/ev2j/wavelength_min+shiftAME

    ############# semi theoretical modelling ############ Uses this method if the material is real
    if(realMaterial):
        #print('real absorption')
        print(fileName)
        data = np.loadtxt(fileName) # Here the absorption coefficient file comes in.
        #data = fileName
        lenData=len(data[:,0]) 
        absM=np.ndarray(shape=(lenData,2), dtype=float) # Creating an array with the same dimension as the file
    # first column energy in eV, second column absorption coefficient in nm**(-1)
        absM[:,0] = data[:,0] # Assigning the energy column in a new array
        absM[:,1] = (data[:,1])/1000. # converting the absorption coefficient from mm⁻1 to nm⁻1 ?
        absM[:,1] =  sigmoid(absM[:,0], Bg, broaden)*absM[:,1] # To smedge the absorption data, to make it nicer.
        
        if(shiftEne >= Bg): # Takes care off the direct bandgap ?
            absM[:,0]=absM[:,0]+(shiftEne-Bg)
            x=np.insert(absM[:,0],0,0) # x is an one dimensional vector, that contains the energy , with 0 as its first object
            y=np.insert(absM[:,1],0,0) # y is an one dimensional vector, that contains the absorption coefficient, with 0 as its first object
        else: # Takes care off the indirect bandgap ?
            absM[:,0]=absM[:,0]+(shiftEne-Bg)
            absM[:,0]=(absM[:,0]>=0)*absM[:,0] # Only the values larger than 0 will remain, so values smaller than 0 becomes 0
            x=absM[:,0] # x is an one dimensional vector, that contains the energy, with 0 as its first object
            y=absM[:,1] # y is an one dimensional vector, that contains the absorption coefficient, with 0 as its first object

        absE_max=x[x.argmax()] # The largest value in x is now assigned to absE_max
    else:
    ############# pure theoretical modelling ############ Uses this method if the material is finctional?
        absE_max=10
        ##### (x-ctPt)**(1./2)/x, 0.001.
        x=(np.linspace(0.001,absE_max,1000001)) # x is an array with values ranging from 0.001 to 10 with 1000001 points in between
        absM=numpy.zeros(shape=(1000001,2)) # creates a two dimensional array with zeros and 1000001 points

        absM[:,0]=x # assigns the first row in absM as x

        

        #######  (x-ctPt)**(1./2)/x, ctPt=0.01, 
        ######  CHANGE 0.000
        ctPt=0.01
        absM[:,1]=5E-4*(x>=ctPt)*(abs(x-ctPt))**(1./2)/x # Only x-values larger or equal than ctPt goes through the calculations, the others becomes 0
        absM[:,0]=absM[:,0]+shiftEne-ctPt

       
        x=np.insert(absM[:,0],0,0) # x is an one dimensional vector, that contains the energy, with 0 as its first object
        y=np.insert(absM[:,1],0,0) # y is an one dimensional vector, that contains the absorption coefficient, with 0 as its first object



    if(relativeRatio):
        #print("relativeRatio is true")
        absE_max=relativeValue
    #else:
        #print("relativeRatio is false")

    ABSinterp = scipy.interpolate.interp1d(x, y) # Interpolating the absorption data with the energy to make the function a(E).

    ###################### total income light  ##############################

# More convenient for to change the units for the solar spectrum, so that one can easily do integrals over photon energy, rather than wavelength, and calculate the number of photons instead of their energy. Therefore, I’ll define the function "SPhotonsPerTEA" which stands for Solar Photons per unit Time, per unit photon Energy-range, per unit Area of the solar cell 

    def SPhotonsPerTEA(E):
        wavelength = plkC*speedL/ev2j/E # E in eV, wavelength in nm
        return AM15interp(wavelength)*(1./E)*(1./E**2) 

    PowerPerTEA = lambda E : E * SPhotonsPerTEA(E) 
    # solar_constant unit: W/m**2
    solar_constant =plkC*speedL*j2ev*scipy.integrate.quad(PowerPerTEA,E_min,E_max,epsabs=1e-100,epsrel=1e-100, full_output=1)[0] # Solar constant" is the sun's total irradiance. If I did this right, it should be 1000 watts/meter 2
	# quad() is ordinary integration; full_output=1 is (surprisingly) how you hide
	# the messages warning about poor accuracy in integrating.


    ############# functions for simpson integral
    def Get_N(a,b,width): # Generates total steps "N" from a start point "a" to a finish point "b" 
    # width is step interval
        N=int((b-a)/width + 1)
        if N%2 == 0:
            N=N+1
        return N

    def GenerateData1(a,b,n,width): # Generates "n" data points for a function from point "a" to point "b"
        datas = []
        r=a
        for i in range(0,n):
            datas.append(func1(r))
            r = r+width
        return datas

    def GenerateData2(a,b,n,width): # Generates "n" data points for a function from point "a" to point "b"
        datas = []
        r=a
        for i in range(0,n):
            datas.append(func2(r))
            r = r+width
        return datas

    def simpson_integral(datas,width,n): # Calculates the integral of a function with a set of data points using the simspon integral 
        sum = datas[0]+datas[n-1]
        for i in range(2,n):
            if i%2== 0:
                sum = sum +4*datas[i-1]
            else:
                sum = sum +2*datas[i-1]
        return sum*width/3.0
    #################### integral cut off #########
    if(relativeRatio): # If there is a relativeratio inserted in the code, this will run, if not this part will be skipped and the "normal" part of the code will run instead.
        ######################  photons above bandgap ######
        def func1(E): 
            wavelength = plkC*speedL/ev2j/E
            return (1-exp(-2*tkCell*ABSinterp(E)))*AM15interp(wavelength)*(1./E)* (1./E**2) # A(E)f_am(E)/E³
        ######################  radiative recombination ####
        def func2(E): 
    #         return ((1-exp(-2.*tkCell*ABSinterp(E)))*E**2./(exp(E/(blmC*tCell*j2ev))-1))
            return ((1-bigfloat.exp(-2.*tkCell*ABSinterp(E),bigfloat.precision(100)))*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(100))-1)) # (A(E)E²)/(exp(ev/kbTc)-1) vs E²/(exp(ev/kbTc)-1)

        a=E_min # dn 
        b=relativeValue

        width=stepE #step
        N=Get_N(a,b,width)
        datas1 = GenerateData1(a,b,N,width) # Generates "N" data points for the function func1 from point "a" to point "b" with step length width
        prefactor0=plkC*speedL*j2ev*qChar*j2ev # h*c*q
        solar_constant1=(prefactor0*simpson_integral(datas1,width,N)) # The "Solar constant" is the sun's total irradiance. If I did this right, it should be 1000 watts/meter2, because that's how NREL normalized their data.

        datas2 = GenerateData2(a,b,N,width) # Generates "N" data points for the function func2 from point "a" to point "b" with step length width
        prefactor1=qChar*2*pi/((speedL/1E9)**2*(plkC*j2ev)**3) # q*2*pi/h³c²
        RR0=(prefactor1*simpson_integral(datas2,width,N)) # Radiative Recombination rate at 0 Quasi-Fermi Level splitting

    #################### normal #########
    else:
        ######################  photons above bandgap ######
        def func(E): 
            wavelength = plkC*speedL/ev2j/E
            return (1-exp(-2*tkCell*ABSinterp(E)))*AM15interp(wavelength)*(1./E)* (1./E**2) # A(E)f_am(E)/E³

        func1=lambda E : func(E) 
        prefactor0=plkC*speedL*j2ev*qChar*j2ev # h*c*q, why is this the prefactor?
        solar_constant1 = prefactor0*scipy.integrate.quad(func1,E_min,E_max,epsabs=1e-100,epsrel=1e-100,full_output=1)[0] # This is the short circuit current density aka Jsc. One can change the integral steps by changing epsrel. 
        #print('')
        #print('-------solar_constant1 A/m^2')
        #print(solar_constant1)

        ######################  radiative recombination ####
        def SPhotonsPerTEA3(E):
            return ((1-bigfloat.exp(-2.*tkCell*ABSinterp(E),bigfloat.precision(500)))*E**2./(bigfloat.exp(E/(blmC*tCell*j2ev),bigfloat.precision(500))-1)) # A(E)E²/exp(E/kbTc) -1
    #        return prefactor1*(1-exp(-2*ABSinterp(E)*tkCell))*E**2/(exp(E/(blmC*tCell))-1)

        integrand = lambda E : SPhotonsPerTEA3(E)
        prefactor1=qChar*2*pi/((speedL/1E9)**2*(plkC*j2ev)**3) # q*2*pi/h³c²
        RR0 = prefactor1*scipy.integrate.quad(integrand, E_min, E_max, epsabs=1e-100,epsrel=1e-100,full_output=1)[0] # J⁰_rad * q. Radiative Recombination rate at 0 QFL splitting. One can change the integral steps by changing epsrel. 
        na=2.51E19 # 1/cm^3
        Eg=gapFun # eV CuSbSeTe
        Cauger=1E-31 # cm^6/s, is the total ambipolar Auger coefficient 
        Ndefect=1E16 # Ndefect = Na, whcich is the acceptor consentration?
        prefactor2=qChar*Cauger*na*na*tkCell*1E-7*1E4*Ndefect*exp(-Eg*ev2j/(blmC*tCell)) # J⁰_aug = q*Ca*d*na²*exp(-Eg/kbTc)*Na, A/m^2, the constants in prefactor 2 & 3 might be wrong or something is written wrong?
        prefactor3=qChar*Cauger*na*na*na*tkCell*1E-7*1E4*exp(-3*Eg*ev2j/(2*blmC*tCell)) # A/m^2, q*Ca*d*na³*exp(-3Eg/2kbTc)*Na ?
    ###################### find maximum  ##############################

# Given what we’ve already done, it’s now simple to calculate the ideal bandgap and efficiency, by numerically maximizing the product JV for each bandgap. 
#The "maximum power point" (MPP) is the point on the JV curve at which this maximum occurs, the maximum power is the power generated at the MPP, and the efficiency is the power divided by the solar constant (i.e. incoming light power).

    def current_density(V):
        # original
        #return (solar_constant1 - RR0*(exp(V*qChar/(blmC*tCell))))
        # paper == clas 1
        #return (solar_constant1 - RR0*(exp(V*qChar/(blmC*tCell))-1)-prefactor2*exp(3*qChar*V/(2*blmC*tCell)))
        return (solar_constant1 - (RR0+1*prefactor2+0*prefactor3*exp(qChar*V/(2*blmC*tCell)))*exp(V*qChar/(blmC*tCell))) # Why is it multiplicated with 0? 
        # clas1

    from scipy.optimize import fmin 

    def fmax(func_to_maximize, initial_guess=0):
	# return the x that maximizes func_to_maximize(x)
        func_to_minimize = lambda x : -func_to_maximize(x)
        return fmin(func_to_minimize, initial_guess, disp=False)[0]

    def V_mpp():
	# voltage at max power point 
        return fmax(lambda V : V * current_density(V))
# Havent maximized the current?

# Calculating the maximum power of the solar cell, by multypling the max voltage and max current density.
    def max_power():
        V = V_mpp()
        return V * current_density(V)

# Calculating the efficiency of the solar cell, by dividing the effect out of the solar cell on the effect from the sun irradiance. 
    def max_efficiency():
        return max_power() / solar_constant

    maxEff=max_efficiency()
    return(maxEff)
    #print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
########### band-gap array #############
    #tkArr=np.array([25,50,100,200,300,400,500,600,700,800,900,1000])
    # band gap of material eV
    BgCISe=0.848945
    BgMat = BgCISe
    # thickness nm
    tkArr=np.array([30, 100, 500, 2000])
    #tkArr=np.linspace(100, 2000, num=11)

    #fileName=['K12.txt','K16.txt','K18.txt','K20.txt','K24.txt','K30.txt']
    # better smooth out the value and starts from exact bandgap engergy
    fileName=['CISe_ABS.txt']

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
        #effDict[fileName[j]]=effTkArr 

    #dFile=pd.DataFrame(effDict,index=tkArr)
    #dFile.to_csv('CIGSe_MAXEFF.txt', header=True, index=True, sep=str(u'\t'), mode='w')
