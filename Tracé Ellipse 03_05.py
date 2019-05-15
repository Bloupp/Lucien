## Dossier de test pour afficher les ellipses en 3D

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
    x= rd*np.cos(phi)*np.cos(teta+v)
    y = rd*np.cos(phi)*np.sin(teta+v)
    z = rd*np.sin(phi)
    return (x,y,z)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#Ya plus qu'a plot le tout en 3 dimenstions
V = [2*pi*i/100 for i in range(100)]
XYZ = [cart(v) for v in V]
X = [x[0] for x in XYZ]
Y = [y[1] for y in XYZ]
Z = [z[2] for z in XYZ]

ax.scatter(X,Y,Z, c='r')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()

#Une fois cela fait, regroupons cela dans un programme

def afficheEllipse(a,e,phi,teta):
    p = a*(1-e**2)
    def r(v):
        return p/(1-e*np.cos(v-teta))
    def cart(v):
        rd = r(v)
        x = rd*np.cos(phi)*np.cos(teta+v)
        y = rd*np.cos(phi)*np.sin(teta+v)
        z = rd*np.sin(phi)
        return [x,y,z]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    V= [pi*i/180 for i in range(360)]
    XYZ = [cart(v) for v in V]
    X = [x[0] for x in XYZ]
    Y = [y[1] for y in XYZ]
    Z = [z[2] for z in XYZ]
    ax.scatter(X,Y,Z, c='r')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()

def affichePlanete(i):
    afficheEllipse(a[i],e[i],incl[i],th[i])

#Enfin réalisons un fonction qui fait appel aux tableaux construits dans la partie sur Le Verrier afin de tracer l'orbite d'une planète à un temps donné

def affichageTrajectoireDuree(i,n,pas):
    phi = incl[i]
    the = th[i]
    dga = a[i]
    def r(v,e):
        return dga*(1-e**2)/(1-e*np.cos(v-the))
    def cart(v,e):
        rd = r(v,e)
        x = rd*np.cos(phi)*np.cos(the+v)
        y = rd*np.cos(phi)*np.sin(the+v)
        z = rd*np.sin(phi)
        return [x,y,z]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    V= [pi*i/180 for i in range(360)]
    for j in range(n):
        exc = excentricite(i,j*pas*annee)
        XYZ = [cart(v,exc) for v in V]
        X = [x[0] for x in XYZ]
        Y = [y[1] for y in XYZ]
        Z = [z[2] for z in XYZ]
        ax.scatter(X,Y,Z)
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()

def affichageSystemeSolaire(t):
    phi,the,dga = 0,0,0
    def r(v,e):
        return dga*(1-e**2)/(1-e*np.cos(v-the))
    def cart(v,e):
        rd = r(v,e)
        x = rd*np.cos(phi)*np.cos(the+v)
        y = rd*np.cos(phi)*np.sin(the+v)
        z = rd*np.sin(phi+v)
        return [x,y,z]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    V= [pi*i/180 for i in range(360)]
    ax.scatter([0],[0],[0])
    for i in range(8):
        phi = incl[i]
        the = th[i]
        dga = a[i]
        exc = excentricite(i,t)
        XYZ = [cart(v,exc) for v in V]
        X = [x[0] for x in XYZ]
        Y = [y[1] for y in XYZ]
        Z = [z[2] for z in XYZ]
        ax.scatter(X,Y,Z)
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()
