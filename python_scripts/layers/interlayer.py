import os
import numpy as np
import matplotlib.pyplot as plt




def lldistance(thisdir):
    
    #filenames = ["DFT_D3.txt", "LDA.txt", "PBE.txt", "rev_vdW_DF2.txt", "vdW_opt88.txt"]
    
    filename = []
    
    # r=root, d=directories, f = files
    for r, d, f in os.walk(thisdir):
        for file in f:
            if file.endswith(".txt"):
                filename.append(os.path.join(r, file))
    
    #print(filename)
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


        temp = etot[1]
        etot.pop(1)
        etot.append(temp)
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



    return etot_0, etot_1, etot_2, etot_3, etot_4


def plotlldistance(etot_0, etot_1, etot_2, etot_3, etot_4):
    
    x = np.linspace(0.0, 5.0, 11)
    n = 24

    plt.plot(x, etot_0/n, "-or")
    plt.plot(x, etot_1/n, "-ob")
    plt.plot(x, etot_2/n, "-og")
    plt.plot(x, etot_3/n, "-oy")
    plt.plot(x, etot_4/n, "-o")
    plt.legend(["PBE", "DFT_D3", "LDA", "rev_vdw", "vdw_opt88"])
    plt.title("All functionals")
    plt.xlabel("Interlayer space [Å]")
    plt.ylabel("Energy/atom [meV]")
    plt.show()

    plt.plot(x, etot_0/n, "-or")
    plt.plot(x, etot_1/n, "-og")
    plt.plot(x, etot_3/n, "-oy")
    plt.plot(x, etot_4/n, "-o")
    plt.legend(["PBE", "DFT_D3", "rev_vdw", "vdw_opt88.txt"])
    plt.title("All functionals except LDA")
    plt.xlabel("Interlayer space [Å]")
    plt.ylabel("Energy/atom [meV]")
    plt.show()


    plt.plot(x, etot_2/n, "-ob")
    plt.legend(["LDA"])
    plt.title("LDA functional")
    plt.xlabel("Interlayer space [Å]")
    plt.ylabel("Energy/atom [meV]")
    plt.show()

