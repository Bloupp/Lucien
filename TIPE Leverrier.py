import numpy as np
from numpy.polynomial import Polynomial
import numpy.linalg as alg
import matplotlib.pyplot as plt
import math
import mpl_toolkits.mplot3d as Axes3D


## Valeurs
G = 6.67E-11
MS = 1.9884E30
MT = 5.972E24
PT = 1.741E11
VP = 30.4E3

#Faire une BDD contenant les valeurs utiles au calcul de (p,v) et [p,v]


## Détermination de A :

#Calcul des coefficient (p,v)
def coef_p(p,v):
    
#Calcul des coefficient [p,v]
def coef_v(p,v):
    
#Calcul de A, l=liste des numéros de planète prises en compte
def A(l):
    n = len(l)
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            
            
#Calcul de A², l=liste des numéros de planète prises en compte
def B(l):
    return np.dot(A,A)
    

##Méthode de détermination des valeurs propres par l'algorithme de Leverrier :

#On fait tout d'abord une fonction permettant de calculer la suite des valeurs Mk que l'on range dans une liste
def M(A):
    n = len(A)
    I = np.eye(n)
    L = [A]
    for i in range(1,n):
        D = L[-1]
        M = np.dot(A,(D-(np.trace(D)/i)*I))
        L.append(M)
    return L

#Ensuite, on calcul le polynome caractéristique avec la méthode annoncée
def Leverrier(A):
    n = len(A)
    L = M(A) #dans cette liste, Mk est à la kieme position
    P = []
    for k in range(n):
        P.append(-np.trace(L[n-k-1])/(n-k))
    P.append(1)
    return Polynomial(P)
    
A=[[1,0,0],[1,1,0],[1,2,1]]

#On va donc calculer les racines de ce polynome qui sont les valeurs propres de A
def valeurs_propres(A):
    p = Leverrier(A)
    return p.roots()


## Résolution de l'équation différentielle

#Diagonalisation : matrice diagonale
def D(l):
    n = len(l)
    D = np.zeros((n,n))
    vp = valeurs_propres( B(l) )
    for i in range(n):
        D[i,i] = vp[i]
    return D
    
#Diagonalisation : P
def P(l):
    
#Définition de K = PH

#Résolution de H :
def K(t)
    

##Tentative de résolution du système d'équations différentielles par la méthode analytique

#La force gravitationnelle est 
def Fg(M,m,r):
    return -G*M*m/(r*r)

'''On va mettre la position de la planète dans un tableau de tableaux contenant les coordonnées cartésiennes au cours du temps, on va se placer dans le référentiel héliocentrique
De plus, on va aussi gérer un tableau contenant la vitesse de l'astre à un instant t'''

#On va calculer la distance entre 2 points
def distance(L,P):
    return np.sqrt((L[0]-P[0])**2+(L[1]-P[1])**2+(L[2]-P[2])**2)

'''Pour une seule planète :
m = masse, X = vecteur position, X_p = vecteur vitesse, t = temps de la modelisation, n = nombre de subdivisions'''
def trajectoire(m,X,X_p,t,n):
    P = [X]
    V = [X_p]
    dt = t/n
    for i in range(n):
        x,y,z=P[-1][0],P[-1][1],P[-1][2]
        vx,vy,vz=V[-1][0],V[-1][1],V[-1][2]
        a_r = Fg(MS,m,distance(P[-1],[0,0,0]))/m #Accélération selon er
        #On projette sur la base cartésienne afin d'obtenir 3 équations différentielles d'ordre 2 mélées
        d = np.sqrt(x**2+y**2+z**2)
        ax = a_r*x/d
        ay = a_r*y/d
        az = a_r*z/d
        #Puis la méthode d'Euler
        P.append([x+vx*dt,y+vy*dt,z+vz*dt])
        V.append([vx+ax*dt,vy+ay*dt,vz+az*dt])
    return np.array(P)

P = trajectoire(MT,[PT,0,0],[0,VP,0],63115200,10000)
X = P[:,0]
Y = P[:,1]
plt.axis('equal')
plt.plot(X,Y)
plt.show()
