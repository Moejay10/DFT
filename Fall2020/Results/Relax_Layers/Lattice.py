import numpy as np
import matplotlib.pyplot as plt

num_layers = 15
a_list = []
b_list = []
c_list = []
layer_list = []

layer_list.append(0)
a_list.append(8.9261)
b_list.append(6.7202)
c_list.append(11.5032)

for i in range(1, num_layers + 1):

    filename = "CONTCAR_" + str(i)
    layer_list.append(i)


    with open(filename) as file:
        lines = file.readlines()
        #Skip the first two lines
        for j in range(2,3):
            line = lines[j]
            pieces = line.split()
            a = float(pieces[0])

        for j in range(3,4):
            line = lines[j]
            pieces = line.split()
            b = float(pieces[1])

        for j in range(4,5):
            line = lines[j]
            pieces = line.split()
            c = float(pieces[2])

    a_list.append(a)
    b_list.append(b)
    c_list.append(c)

a_list = np.array(a_list)
b_list = np.array(b_list)
c_list = np.array(c_list)
layer_list = np.array(layer_list)

plt.plot(layer_list, a_list, "-or")
plt.title("Lattice Constant a")
plt.show()

plt.plot(layer_list, b_list, "-og")
plt.title("Lattice Constant b")
plt.show()

plt.plot(layer_list, c_list, "-ob")
plt.title("Lattice Constant c")
plt.show()
