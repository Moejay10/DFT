import numpy as np
import os



def readFilesFromDirectory(thisdir):
    
    
    filepaths = []
    filenames = []
    
    # r=root, d=directories, f = files
    for r, d, f in os.walk(thisdir):
        for file in f:
            if file.endswith(".txt"):
                filepaths.append(os.path.join(r, file))
                filenames.append(file)
    
    return filepaths, filenames

def readEnergyFromFile(filepaths, filenames):
    
    data = {}
    j = 0
    for files in filepaths:
        file = open(files, "r")
        lines_read = file.readlines()
        N = len(lines_read)

        etot = []

        for i in range(1,N):
            line = lines_read[i]
            pieces = line.split()
            a = float(pieces[2])
            etot.append(a)

        data[filenames[j]] = np.array(etot)
        j += 1


    return data



