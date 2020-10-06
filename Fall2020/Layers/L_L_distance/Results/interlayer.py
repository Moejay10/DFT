import os
import numpy as np
import matplotlib.pyplot as plt


filename = ["DFT_D3.txt", "LDA.txt", "PBE.txt", "rev_vdw.txt", "vdw_opt88.txt"]

for files in filename:
    file = open(files, "r")
    lines_read = file.readlines()
    N = len(lines_read)

    etot = []

    for i in range(1,N):
        line = lines_read[i]
        pieces = line.split()
        a = float(pieces[2])
        etot.append(a)

    if files == filename[0]:
        etot_0 = np.array(etot)

    elif files == filename[1]:
        etot_1 = np.array(etot)

    elif files == filename[2]:
        etot_2 = np.array(etot)

    elif files == filename[3]:
        etot_3 = np.array(etot)

    else:
        etot_4 = np.array(etot)

x = np.linspace(0.05, 0.05*N, N)
n = 24
print(x)

plt.plot(x[:-1], etot_0/n, "-or")
plt.plot(x[:-1], etot_1/n, "-ob")
plt.plot(x[:-1], etot_2/n, "-og")
plt.plot(x[:-1], etot_3/n, "-oy")
plt.plot(x[:-1], etot_4/n, "-o")
plt.legend(["DFT_D3.txt", "LDA.txt", "PBE.txt", "rev_vdw.txt", "vdw_opt88.txt"])
#plt.title("")
plt.xlabel("Interlayer space [Å]")
plt.ylabel("Energy/atom [meV]")
plt.show()
