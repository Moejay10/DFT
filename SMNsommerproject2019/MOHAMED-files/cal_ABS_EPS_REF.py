#!/software/apps/python/2.7.6/snic-1/bin/python

#definition for the function of the absorption
def absimple(E, EPSR, EPSI):
        from math import sqrt
        # i need E for energy and EPSI and EPSR
        # con= h*c/2*pi
        con=1.973187364683e-5
        ab=sqrt(2.*(E**2)*(-EPSR+sqrt(EPSI*EPSI+EPSR*EPSR ))/(con**2) )/10000.
        return ab

#import numpy
import numpy
from math import *

#the "keylines" and "keyword" gave me the number to start
# and the position to start respectively
keylines= 'NEDOS'
keyword= 'IMAGINARY DIELECTRIC'

f=open('OUTCAR','r')
#with open('datafile.3') as f:
lineno = 0
templine=numpy.array([0])
tempposition=numpy.array([0])
while True:
   line=f.readline()
   lineno = lineno+1
   if not line: break
   if keylines in line:
         fields=line.split('  ')
     a=int(fields[5])
   if keyword in line:
    print repr(keyword), 'found in line ', lineno, 'in datafile :',repr(line)
    templine=numpy.concatenate(([lineno],templine))
    tempposition=numpy.concatenate(([f.tell()],tempposition))
print a
print lineno
print templine 
print tempposition

epsi = numpy.zeros((a,7))
epsr = numpy.zeros((a,7))
EPSfinal = numpy.zeros((a,7))
absarray = numpy.zeros((a,4))
absarray2C = numpy.zeros((a,2))
refractiveArray = numpy.zeros((a-1,3))

   
print(f.tell())
f.seek(tempposition[0],0)
print(f.tell())

counter=0
temp=0
while True:
 line = f.readline()
 counter=counter+1
 if (counter > 2)&( counter < (a+3) ): 
   epsi[counter-3] = numpy.fromstring(line,dtype=float,sep='  ')
 

 if (counter > a+6 ) & (counter< a+a+7):
     epsr[counter-a-7]=numpy.fromstring(line,dtype=float,sep='  ')
 if counter>(a+a +6):
   break

f.close()
    
ab=numpy.vectorize(absimple)

absarray[:,0]= epsi[:,0]
absarray[:,1]= ab(epsi[:,0],epsr[:,1],epsi[:,1])
absarray[:,2]= ab(epsi[:,0],epsr[:,2],epsi[:,2])
absarray[:,3]= ab(epsi[:,0],epsr[:,3],epsi[:,3])

absarray2C[:,0]= epsi[:,0]
absarray2C[:,1]= (absarray[:,1] + absarray[:,2] + absarray[:,3])/3.0

EPSfinal= numpy.column_stack((epsi[:,0],epsr[:,1],epsi[:,1],epsr[:,2],epsi[:,2],epsr[:,3],epsi[:,3]))



tempR=(epsr[1:,1]+epsr[1:,2]+epsr[1:,3])/3.0
tempI=(epsi[1:,1]+epsi[1:,2]+epsi[1:,3])/3.0

tempN=numpy.sqrt((numpy.sqrt(tempR**2+tempI**2)+tempR)/2.0)
tempK=numpy.sqrt((numpy.sqrt(tempR**2+tempI**2)-tempR)/2.0)
waveLength=1240/epsi[1:,0]
refractiveArray = numpy.column_stack((waveLength, tempN, tempK))

numpy.savetxt('EPS.txt',(EPSfinal),'%10.6f   %10.6f   %10.6f   %10.6f   %10.6f   %10.6f   %10.6f','\n' )

numpy.savetxt('ABS.txt',(absarray),'%10.6f   %10.6f   %10.6f   %10.6f','\n')

numpy.savetxt('ABS_oneCol.txt',(absarray2C),'%10.6f   %10.6f','\n')

numpy.savetxt('Refractive.txt',(refractiveArray),'%10.6f   %10.6f   %10.6f','\n')




