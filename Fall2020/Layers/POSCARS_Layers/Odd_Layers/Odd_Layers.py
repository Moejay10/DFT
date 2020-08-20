import numpy as np
import os
import fileinput


filename = input("Write the filename here: ")
num_lines_to_delete = int(input("How many lines do you want to delete: "))


a = 0
with open(filename) as file:
    lines = file.readlines()
    #Skip the first two lines
    for j in range(2,3):
        line = lines[j]
        pieces = line.split()
        b = pieces[0]
        a = float(pieces[0]) - 4.504734039

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

with open(filename,"r") as f:
    newline=[]
    for word in f.readlines():
        newline.append(word.replace(ba,Ba))


with open(filename,"w") as f:
    for line in newline:
        f.writelines(line)

with open(filename,"r") as f:
    newline=[]
    for word in f.readlines():
        newline.append(word.replace(si,Si))


with open(filename,"w") as f:
    for line in newline:
        f.writelines(line)

for i in range(num_lines_to_delete):
    line_to_delete = int(input("Type the line you want to delete: ")) - i + 8
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



new_filename = input("Write the new filename here: ")

os.rename(filename, new_filename)
