import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


#---------------------------------1Layer---------------------------------------#


ENCUT = np.array([300, 350, 400, 450, 500, 550, 600, 650])

Etot_Layer1 = np.array([-32.508317, -32.521323, -32.522743, -32.5228, -32.525675,
-32.523672,
-32.523697,
-32.527122])

sns.set()
plt.title("Convergence of Etot for 1Layer")
plt.plot(ENCUT, Etot_Layer1)
plt.xlabel("ENCUT [eV]")
plt.ylabel("Etot [eV]")
plt.show()

E_relative_1Layer = np.array([0.000157000000002,
3.20000000044729E-05,
0.000230999999999,
0.000295999999999,
0.000588999999998,
0.000363999999998,
0.000714000000002
])

sns.set()
plt.title("Convergence of Erelative for 1Layer")
plt.plot(ENCUT[1:], E_relative_1Layer)
plt.xlabel("ENCUT [eV]")
plt.ylabel("E_relative [eV]")
plt.show()


k_density_1Layer = np.array([0.007803000000003,
4.89999999970792E-5,
0.000172999999997,
1.40000000001805E-5,
5.99999999906231E-6,
2.99999999953116E-6,
0
])

k_points = np.array([2,
3,
4,
5,
6,
7,
8,
9
])

sns.set()
plt.title("Convergence of k-point density for 1Layer")
plt.plot(k_points[:-1], k_density_1Layer)
plt.xlabel("k-points []")
plt.ylabel("E_relative [eV]")
plt.show()

Vacuum_Energy_1Layer = np.array([-30.842645,
-34.720515,
-33.149405,
-32.697511,
-32.576047,
-32.542593,
-32.531414,
-32.527677,
-32.524935,
-32.523865,
-32.523892
])

Vacuum_Space_Layer = np.array([0,
2,
4,
6,
8,
10,
12,
14,
16,
18,
20
])


sns.set()
plt.title("Convergence of Vacuum for 1Layer")
plt.plot(Vacuum_Space_Layer, Vacuum_Energy_1Layer)
plt.xlabel("Vacuum [Å]")
plt.ylabel("Etot [eV]")
plt.show()

#---------------------------------2Layers--------------------------------------#

Etot_2Layers = np.array([-72.382356,
-72.408616,
-72.410879,
-72.412681,
-72.418037,
-72.413211,
-72.41364,
-72.417776
])

sns.set()
plt.title("Convergence of Etot for 2Layers")
plt.plot(ENCUT, Etot_2Layers)
plt.xlabel("ENCUT [eV]")
plt.ylabel("Etot [eV]")
plt.show()

E_relative_2Layers = np.array([0.00038099999999,
0.000352000000007,
0.00048799999999,
0.001701999999995,
0.000918999999996,
0.000981999999993,
0.000726
])

sns.set()
plt.title("Convergence of Erelative for 2Layers")
plt.plot(ENCUT[1:], E_relative_2Layers)
plt.xlabel("ENCUT [eV]")
plt.ylabel("E_relative [eV]")
plt.show()


k_density_2Layers = np.array([0.000461999999999,
4.80000000067093E-05,
3.80000000035352E-05,
1.60000000022364E-05,
1.99999999495049E-06,
2.99999999242573E-06
])



sns.set()
plt.title("Convergence of k-point density for 2Layers")
plt.plot(k_points[1:-1], k_density_2Layers)
plt.xlabel("k-points []")
plt.ylabel("E_relative [eV]")
plt.show()

Vacuum_Energy_2Layers = np.array([-81.126685,
-75.326125,
-73.155085,
-72.622354,
-72.481231,
-72.43981,
-72.425718,
-72.420102,
-72.416867,
-72.414742,
-72.413832
])




sns.set()
plt.title("Convergence of Vacuum for 2Layers")
plt.plot(Vacuum_Space_Layer, Vacuum_Energy_2Layers)
plt.xlabel("Vacuum [Å]")
plt.ylabel("Etot [eV]")
plt.show()
