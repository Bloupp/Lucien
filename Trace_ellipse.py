import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

pi = np.pi

e = 0.20564
a = 57909227000
phi = pi*7/180
teta = pi*48/180

#Considérons une ellipse inclinée d'un angle teta sur le plan equatorial et d'un angle phi sur son orthogonal
# On définit tout d'abord p
p = a*(1-e**2)

#Ensuite on calcule r pour un angle sur l'ellipse donné v
def r(v):
    return p/(1-e*np.cos(v-teta))

# On calcul ensuite les coordonnées cartésiennes :
def cart(v):
    rd = r(v)
    x= rd*np.cos(phi)*np.cos(teta)
    y = rd*np.cos(phi)*np.sin(teta)
    z = rd*np.sin(phi)
    return (x,y,z)

#Ya plus qu'a plot le tout en 3 dimenstions
V = [2*pi*i/100 for i in range(100)]
X = [cart(v)[0] for v in V]
Y = [cart(v)[1] for v in V]
Z = [cart(v)[2] for v in V]

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(X,Y,V)
ax.legend()

plt.show()
