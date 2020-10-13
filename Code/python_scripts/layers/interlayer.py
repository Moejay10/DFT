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
            pbe = np.array(etot)

        elif files == filename[1]:
            dft_d3 = np.array(etot)

        elif files == filename[2]:
            lda = np.array(etot)

        elif files == filename[3]:
            rev = np.array(etot)

        else:
            vdw = np.array(etot)



    return pbe, dft_d3, lda, rev, vdw



