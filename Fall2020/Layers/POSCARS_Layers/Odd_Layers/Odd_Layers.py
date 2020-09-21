import numpy as np
import os
import fileinput

"""
NOTE: POSCAR-4L does not work
"""



layer = input("Write which layer here: ")
filename = "POSCAR_"  + layer + "L"
num_lines_to_delete = 12


N_list1 = filename.strip().split("_")
N_list2 = N_list1[1].strip().split("L")
N = int(N_list2[0])

list = []

for i in range(1, num_lines_to_delete + 1):
    if i == 7 or i == 8:
        list.append(N*i - int(N/2))
    else:
        list.append(N*i)

print(list)


a = 0
with open(filename) as file:
    lines = file.readlines()
    #Skip the first two lines
    for j in range(2,3):
        line = lines[j]
        pieces = line.split()
        b = pieces[0]
        a = float(pieces[0]) - 4.5047340393

a = float("{:.10f}".format(a))
a = str(a)

with open(filename,"r") as f:
    newline=[]
    for word in f.readlines():
        newline.append(word.replace(b,a))


with open(filename,"w") as f:
    for line in newline:
        f.writelines(line)

with open(filename) as file:
    lines = file.readlines()
    #Skip the first two lines
    for j in range(6,7):
        line = lines[j]
        pieces = line.split()
        ba = pieces[0]
        si = pieces[1]
        Ba = int(pieces[0]) - 4
        Si = int(pieces[1]) - 8

Ba = str(Ba)
Si = str(Si)

the_file = open(filename, "r")
lines_read = the_file.readlines()
lines_read[6] = Ba + "     " + Si + "\n"

the_file = open(filename, "w")
the_file.writelines(lines_read)
the_file.close()




i = 0
while i < num_lines_to_delete:
    line_to_delete = list[i] - i + 8
    initial_line = 1
    file_lines = {}

    with open(filename) as f:
        content = f.readlines()

    for line in content:
        file_lines[initial_line] = line.strip()
        initial_line += 1

    f = open(filename, "w")
    for line_number, line_content in file_lines.items():
        if line_number != line_to_delete:
            f.write('{}\n'.format(line_content))

    f.close()
    print('Deleted line: {}'.format(line_to_delete))

    i += 1

new_layer = int(layer) - 1
new_filename = "POSCAR_" + "L" + str(new_layer)

os.rename(filename, new_filename)
