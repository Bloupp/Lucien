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
#On va avoir un tableau m,n,a,N,P contenant les valeurs

## Détermination de A :

#Calcul des coefficient (p,v)
def coef_p(p,v):
    return (2*f*m[v]*N[p,v]/(n[p]*(a[p]**2))
    
#Calcul des coefficient [p,v]
def coef_c(p,v):
    return (2*f*m[v]*P[p,v]/(n[p]*(a[p]**2))
 
#Construction de A, l=liste des numéros de planètes prises en compte
def A(l):
    n = len(l)
    A = np.zeros((n,n))
    for i in l:
        for j in l:
            if i=j then:
                s_i = 0
                for k in l:
                    if k != i then:
                        s += coef_p(i,k)
                A[i,i] = s
            else:
                A[i,j] = coef_c(i,j)
    return A
                
#Calcul de A², l=liste des numéros de planètes prises en compte
def B(l):
    return np.dot(A(l),A(l))
    

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
def K(t):
    
            
            
