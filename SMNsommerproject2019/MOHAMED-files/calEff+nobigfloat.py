from __future__ import division, print_function, unicode_literals
import numpy, scipy.interpolate, scipy.integrate
from pylab import *
import matplotlib.pyplot as plt
#import bigfloat
import time

import pandas as pd

import sys

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
    dataAM=np.loadtxt('/Users/rongzhenchen/workplace/bookChapter/DIELECTRICFUNCTION/pd/data/AM.txt')

    # band-gap of loaded material (eV)

    #kb=0.02585...
    broaden=0.000001
    #    broaden=0.01
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

    ############# semi theoretical modelling ############
    if(realMaterial):
        #print('real absorption')
        print(fileName)
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
        '''
        absM[:,1]=5E-4*(x>=shiftEne)
        absM[:,0]=absM[:,0]
        '''

        #######  (x-ctPt)**(1./2)/x, ctPt=0.01, 
        ######  CHANGE 0.000
        ctPt=0.01
        absM[:,1]=5E-4*(x>=ctPt)*(abs(x-ctPt))**(1./2)/x
        absM[:,0]=absM[:,0]+shiftEne-ctPt

        #######  (x)**(2)
        '''
        absM[:,1]=5E-4*(x)**(4./2)
        absM[:,0]=absM[:,0]+shiftEne
        '''

        #######  exp(-x)
        '''
        ctPt=0.2
        absM[:,1]=25E-4*(x>=ctPt)*exp(-x)
        absM[:,0]=absM[:,0]+shiftEne-ctPt
        '''

        x=np.insert(absM[:,0],0,0)
        y=np.insert(absM[:,1],0,0)


     ######### test 
     ######### f=C*exp(-(E-Ec).^2/Ec)   
     ######### a(E)=a(E).*f


    '''
    C=20

    Ec=Bg+0.5
    fE=C*exp(-(x-Ec)**2/Ec)
    new=fE*y
    plot(x,y,label='orig ABS')
    y=new
    plot(x,new,label="new ABS")
    plt.legend( )
    plt.show()
    plot(x,y,label='orig ABS')
    arrEEE=np.array([-0.5,0,0.5,1])
    for i in range(0,4):
       Ec=Bg+arrEEE[i]
       strLeg='Bg+'+str(arrEEE[i])
       fE=C*exp(-(x-Ec)**2/Ec)
       new=fE*y
       axarr[1].plot(x,new,label=strLeg,linewidth=2.0)

    Ec=2.50
    strLeg='peak of AM1.5'
    fE=C*exp(-(x-Ec)**2/Ec)
    new=fE*y
    axarr[1].plot(x,new,label=strLeg,linewidth=2.0)
    axarr[1].legend(loc=3)
    axarr[1].set_xlabel(' Energy [eV]')
    axarr[1].set_ylabel(' ABS [10^4 cm^-1]')
    plt.show()
    '''

    if(relativeRatio):
        #print("relativeRatio is true")
        absE_max=relativeValue
    #else:
        #print("relativeRatio is false")

    ABSinterp = scipy.interpolate.interp1d(x, y)

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
            return (1-exp(-2*tkCell*ABSinterp(E)))*AM15interp(wavelength)*(1./E)* (1./E**2)
        ######################  radiative recombination ####
        def func2(E): 
    #         return ((1-exp(-2.*tkCell*ABSinterp(E)))*E**2./(exp(E/(blmC*tCell*j2ev))-1))
            return ((1-exp(-2.*tkCell*ABSinterp(E)))*E**2./exp(E/(blmC*tCell*j2ev))-1))

        a=E_min # dn 
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
            return (1-exp(-2*tkCell*ABSinterp(E)))*AM15interp(wavelength)*(1./E)* (1./E**2)

        func1=lambda E : func(E) 
        prefactor0=plkC*speedL*j2ev*qChar*j2ev
        solar_constant1 = prefactor0*scipy.integrate.quad(func1,E_min,E_max,epsabs=1e-100,epsrel=1e-100,full_output=1)[0] 
        #print('')
        #print('-------solar_constant1 A/m^2')
        #print(solar_constant1)

        ######################  radiative recombination ####
        def SPhotonsPerTEA3(E):
            return ((1-exp(-2.*tkCell*ABSinterp(E)))*E**2./(exp(E/(blmC*tCell*j2ev))-1))
    #        return prefactor1*(1-exp(-2*ABSinterp(E)*tkCell))*E**2/(exp(E/(blmC*tCell))-1)

        integrand = lambda E : SPhotonsPerTEA3(E)
        prefactor1=qChar*2*pi/((speedL/1E9)**2*(plkC*j2ev)**3)
        RR0 = prefactor1*scipy.integrate.quad(integrand, E_min, E_max, epsabs=1e-100,epsrel=1e-100,full_output=1)[0]
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
        return (solar_constant1 - (RR0+1*prefactor2+0*prefactor3*exp(qChar*V/(2*blmC*tCell)))*exp(V*qChar/(blmC*tCell)))
        # clas1

    from scipy.optimize import fmin 

    def fmax(func_to_maximize, initial_guess=0):
        func_to_minimize = lambda x : -func_to_maximize(x)
        return fmin(func_to_minimize, initial_guess, disp=False)[0]

    def V_mpp():
        return fmax(lambda V : V * current_density(V))

    def max_power():
        V = V_mpp()
        return V * current_density(V)

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
