import numpy as np
import os
import fileinput
from shutil import copyfile
from sys import exit
import argparse


def main():
    
    parser = argparse.ArgumentParser(
    description='Constructs interlayer structure with different interlayer distance.',
    formatter_class=argparse.RawTextHelpFormatter
    )


    parser.add_argument("-i", "--input", type=str, help="The folder name of where the initial structures is located.", choices=['BaSi2', 'Graphite'], required=True)

    parser.add_argument("-o", "--output", type=str, help="The folder name of where the output files are going to be stored.", required=True)

    parser.add_argument("-s", "--scale", type=float, help="Scale factor to decide interlayer distances", required=True)


    args = parser.parse_args()

    src = str(walk_up_folder(os.getcwd(), 4)) + "/Data/Layer/L_L_distance/" + args.input + "/POSCARS/" 
    
    if args.input == 'BaSi2':
        construct_POSCAR(src, basi2, args.output, args.scale)

    elif args.input == 'Graphite':
        construct_POSCAR(src, graphite, args.output, args.scale)
    
    else:
        print("Wrong command line argument. Input must either be 'BaSi2' or 'Graphite'.")

def construct_POSCAR(src, func, output, scale):
        
        folders = ["/DFT-D3/", "/PBE/", "/LDA/", "/rev-vdW-DF2/", "/vdW-opt88/"]
        dst = src + "outputFiles/" + output 
        src = src + "inputFiles/"
        
        for folder in folders:
            src1 = src + folder
            dst1 = dst + folder
            
            if not os.path.isdir(dst1):
                os.makedirs(dst1)

            func(src1, dst1, scale)


def walk_up_folder(path, depth=1):
    _cur_depth = 1
    while _cur_depth < depth:
        path = os.path.dirname(path)
        _cur_depth += 1
    return path


def basi2(src, dst, scale):
         
    layer = 2
    old_file = src + "CONTCAR"
    num_lines_to_change = 12

    N = 2

    list = []

    for i in range(1, num_lines_to_change + 1):
        if i == 7 or i == 8:
            list.append(N*i - int(N/2))
        else:
            list.append(N*i)

    num_files = 10

    for i in range(1, num_files + 1):


        new_layer = str(i)
        new_file = dst +  "/POSCAR_" + new_layer

        source = old_file
        target = new_file

        # adding exception handling
        try:
            copyfile(source, target)
        except IOError as e:
            print("Unable to copy file. %s" % e)
            exit(1)
        except:
            print("Unexpected error:", sys.exc_info())
            exit(1)

        #print("\nFile copy done!\n")

        """
        while True:
            print("Do you like to print the file ? (y/n): ")
            check = input()
            if check == 'n':
                break
            elif check == 'y':
                file = open(target, "r")
                print("\nHere follows the file content:\n")
                print(file.read())
                file.close()
                print()
                break
            else:
                continue
        """

        
        with open(new_file) as file:
            lines = file.readlines()
            #Skip the first two lines
            for j in range(2,3):
                line = lines[j]
                pieces = line.split()
                old = pieces[0]
                new = float(pieces[0]) + scale*i
            
            new = format(new, '.16f')
        
        with open(new_file,"r") as f:
            newline=[]
            for word in f.readlines():
                newline.append(word.replace(old,new))


        with open(new_file,"w") as f:
            for line in newline:
                f.writelines(line)

        """ 

        the_file = open(new_file, "r")
        lines_read = the_file.readlines()
        N = len(lines_read)
        #Skip the first eight lines
        count = 0
        for j in range(9, N):
            if j - 9 == list[count] - 1:
                line = lines_read[j]
                pieces = line.split()
                a = float(pieces[0])
                b = (pieces[1])
                c = (pieces[2])
                d = (pieces[3])
                e = (pieces[4])
                f = (pieces[5])

                #a += (scale*i)/2
                #a = float("{:.10f}".format(a))

                lines_read[j] = "     " + str(a) + "         " + b + "         " + c + " " + d + " " + e + " " + f + "\n"
                count += 1


        the_file = open(new_file, "w")
        the_file.writelines(lines_read)
        the_file.close()
        """

def graphite(src, dst, scale):

             
    layer = 2
    old_file = src + "CONTCAR"
    num_lines_to_change = 12

    N = 2

    list = [9, 11]

    num_files = 10

    for i in range(1, num_files + 1):


        new_layer = str(i)
        new_file = dst +  "/POSCAR_" + new_layer

        source = old_file
        target = new_file

        # adding exception handling
        try:
            copyfile(source, target)
        except IOError as e:
            print("Unable to copy file. %s" % e)
            exit(1)
        except:
            print("Unexpected error:", sys.exc_info())
            exit(1)

        #print("\nFile copy done!\n")
        
        with open(new_file) as file:
            lines = file.readlines()
            #Skip the first two lines
            for j in range(4,5):
                line = lines[j]
                pieces = line.split()
                old = pieces[2]
                new = float(pieces[2]) + scale*i

            new = format(new, '.16f')
            new = str(new)
        
        
        with open(new_file,"r") as f:
            newline=[]
            for word in f.readlines():
                newline.append(word.replace(old,new))


        with open(new_file,"w") as f:
            for line in newline:
                f.writelines(line)

        
        """
        the_file = open(new_file, "r")
        lines_read = the_file.readlines()
        N = len(lines_read)
        for j in range(N):
            if (j in list):
                line = lines_read[j]
                pieces = line.split()
                a = (pieces[0])
                b = (pieces[1])
                c = float(pieces[2])
                d = (pieces[3])
                e = (pieces[4])
                f = (pieces[5])

                #c += (scale*i)/2
                #c = float("{:.10f}".format(c))

                lines_read[j] = "     " + a + "         " + b + "         " + str(c) + " " + d + " " + e + " " + f + "\n"


        the_file = open(new_file, "w")
        the_file.writelines(lines_read)
        the_file.close()
        """


if __name__ == '__main__':
    main()
