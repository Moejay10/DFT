import numpy as np
import os
import fileinput
from shutil import copyfile
from sys import exit



layer = 2
old_file = "POSCAR_0"
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
    new_file = "POSCAR_" + new_layer

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

    print("\nFile copy done!\n")

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

            a += 0.5*i
            a = float("{:.10f}".format(a))

            lines_read[j] = "     " + str(a) + "         " + b + "         " + c + " " + d + " " + e + " " + f + "\n"
            count += 1


    the_file = open(new_file, "w")
    the_file.writelines(lines_read)
    the_file.close()
