from matplotlib.pyplot import *
from numpy import *


# Cohesive energy
print("\nCohesive energy:\n")

BaSi2 = [-114.566485, -122.507658, -110.810402, -97.448911, -101.371068, -100.145623]
BaGe2 = [-104.17010734, -114.03012907, -100.08365249, -85.58778651, -88.08850118, -87.77319843]
BaSiC =[-131.47041368, -142.39298132, -128.37758313, -116.08448788, -120.83039681, -119.58352109]
SrSi2 =[-112.39986107, -119.42121147, -108.14423978, -93.81267291, -97.46257342, -96.27520013]



Ba = [-0.06590638, -0.04124352, -0.04124352,  1.28100863, 1.13524021, 1.29693127]
Si = [-0.01645681, -0.01545810, -0.01545810, 0.28888477, -0.12731162,  0.07246815]
Ge = [-0.01304211, -0.01186089, -0.01186089, 0.42525342, 0.08006653, 0.26159535]
C = [0.01410711, 0.02128251, 0.02128251, 0.14417636, -0.42315797, -0.12109662]
Sr = [-0.08443346, -0.07118907, -0.07118907, 1.29865145, 1.17331491, 1.32999537]

N = 24 # number of atoms


def Ecohesive(AB2, A, B2):
    for i in range(len(AB2)):

        answ = (AB2[i] -(8*A[i] + 16*B2[i]))/N

        if i ==0:
            print ("DFT-D3: %g eV\n" % answ)

        elif i ==1:
            print ("LDA: %g eV\n" % answ)

        elif i ==2:
            print ("PBE: %g eV\n" % answ)

        elif i ==3:
            print ("rew-vdW-DF2: %g eV\n" % answ)

        elif i ==4:
            print ("vdW-DF-cx: %g eV\n" % answ)

        elif i ==5:
            print ("vdW-opt88: %g eV\n" % answ)


def Ecohesive2(AB2, A, B, C):
    for i in range(len(AB2)):

        answ = (AB2[i] -(8*A[i] + 8*B[i] + 8*C[i]))/N

        if i ==0:
            print ("DFT-D3: %g eV\n" % answ)

        elif i ==1:
            print ("LDA: %g eV\n" % answ)

        elif i ==2:
            print ("PBE: %g eV\n" % answ)

        elif i ==3:
            print ("rew-vdW-DF2: %g eV\n" % answ)

        elif i ==4:
            print ("vdW-DF-cx: %g eV\n" % answ)

        elif i ==5:
            print ("vdW-opt88: %g eV\n" % answ)


print("BaSi2")
Ecohesive(BaSi2, Ba, Si)
print("\n")

print("BaGe2")
Ecohesive(BaGe2, Ba, Ge)
print("\n")

print("BaSiC")
Ecohesive2(BaSiC, Ba, Si, C)
print("\n")

print("SrSi2")
Ecohesive(SrSi2, Sr, Si)
print("\n")



# Entalphy bonding energy

basi2 = [-97.448911, -110.810402]

ba = [-1.82065575, -12.28668927]
si = [-52.30009236, -57.49310413]
si_atom = [0.28888477, -0.01545810]
ba_atom = [1.28100863, -0.04124352]

N_Ba_Si = 24 # Number of atoms
N_Si = 16 # Number of atoms
N_Ba = 8

print("Entalphy bonding energy\n")
def Hb(AB2, A, B2, Nb):
    for i in range(len(AB2)):

        answ = (AB2[i] -(A[i] + B2[i]))/Nb
        if i ==0:
            print ("rew-vdW-DF2: %g eV\n" % answ)

        elif i ==1:
            print ("PBE: %g eV\n" % answ)

Hb(basi2, ba, si, N_Ba_Si)


a = 1
Hb_Si = (si[a] - 16*si_atom[a])/N_Si
print("\n")
print("Hb Si")
print("%g" % Hb_Si)



Hb_Ba = (ba[a] - 8*ba_atom[a])/N_Ba
print("\n")
print("Hb Ba")
print("%g" % Hb_Ba)




# Entahlpy decomposition
print("\nEntahlpy decomposition \n")
BaSi2 = -97.448911
BaGe2 = -85.58778651
BaSiC = -116.08448788
SrSi2 = -93.81267291

Ba_bulk = -1.35741360
Si_bulk = -41.348542
Ge_bulk = -15.482757
C_bulk = -36.239670
Sr_bulk = -1.226760

Hd1 = (BaSi2 - (4*Ba_bulk + 2*Si_bulk))/24.

Hd2 = (BaGe2 - (4*Ba_bulk + 4*Ge_bulk))/24.

Hd3 = (BaSiC - (4*Ba_bulk + Si_bulk + 2*C_bulk))/24.

Hd4 = (SrSi2 - (2*Sr_bulk + 2*Si_bulk))/24.

print("BaSi2 - Hd")
print("%g" % Hd1)
print("\n")

print("BaGe2 - Hd")
print("%g" % Hd2)
print("\n")

print("BaSiC - Hd")
print("%g" % Hd3)
print("\n")

print("SrSi2 - Hd")
print("%g" % Hd4)
print("\n")


# Bulk modulus

BaSi2 = [514.4973, 690.0785, 203.1588, 228.9940, 203.1588, 244.7660]
BaGe2 = [500.2059, 616.9678, 575.6744, 192.0504, 195.1529, 185.3042]
BaSiC = [-1.4626, 190.6109, 566.6616, 51.3093, 200.6984, 268.7953]
SrSi2 = [570.0940, 779.1431, 682.0571, 244.1555, 246.1233, 273.3157]

#B_0 = 1.0/9.0 *((C_11 + C_22 + C_33) + 2*(C_12 + C_23 + C_31))

print ("Bulk Modulus\n")
def Bulk_Modulus(A):
    B_0 = 1.0/9.0 *((A[0] + A[1] + A[2]) + 2*(A[3] + A[4] + A[5]))
    print("%g\n" % B_0)

print("BaSi2")
Bulk_Modulus(BaSi2)

print("BaGe2")
Bulk_Modulus(BaGe2)

print("BaSiC")
Bulk_Modulus(BaSiC)

print("SrSi2")
Bulk_Modulus(SrSi2)
