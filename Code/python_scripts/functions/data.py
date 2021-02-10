import numpy as np
import pandas as pd
import os
from ase import io


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
    
    df = pd.DataFrame({ key:pd.Series(value) for key, value in data.items() })

    return df

def read_table(filepaths, filenames):
    
    dic = {}
    for file, filename in zip(filepaths, filenames):
        temp_df = pd.read_table(file)
        temp_df['Functional'] = temp_df['Total energy (eV)'].str.extract(r'([\w]+-?[\w]+-?[\w]+)\/O')
        df = temp_df[['Functional', 'E0']]
        dic[filename] = df
    

    return dic


def read_atoms(filepaths):

   atoms = io.read(filepaths)
   # Structure data as a 3x3 matrix
   cell = atoms.cell 

   # Lattice constants
   a = cell[0][0]
   b = cell[1][1]
   c = cell[2][2]

   lattice = [a, b, c]

   return lattice


if __name__ == '__main__':
    
    thisdir = '../../../Results/Bulk/BaGe2/Etot/'
    filepaths, filenames = readFilesFromDirectory(thisdir)
    
    tmp_dic = read_table(filepaths, filenames)
    keys = list(tmp_dic.keys())
    print(tmp_dic[keys[0]])

