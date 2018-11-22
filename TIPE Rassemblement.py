import numpy as np
from numpy.polynomial import Polynomial
import numpy.linalg as alg
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D

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


## Diagonalisation de la matrice M

def D(M):
    n = len(M)
    D = np.zeros((n,n))
    vp = valeurs_propres( M )
    for i in range(n):
        D[i,i] = vp[i]
    return D
    
def Q(M):
    n = len(l)
    vp = valeurs_propres( M )
    X0 = np.zeros((1,n))
    Q = np.zeros((n,0))
    for i in range(n):
        S = M-vp[i]*np.eye(n)
        X = np.linalg.solve(S,X0)
        Q = np.concatenate(Q,X,axis=1)
    return Q

def diagonalisation(M):
    Q=Q(M)
    return inv.alg(Q),D(M),Q
    
    
## Je fais une fonction qui va renvoyer la matrice colonne des fonctions de H et L

#La fonction qui suit se place comme un moyen d'application de l'algorithme construit jusqu'ici, celle-ci est faite dans le but de s'executer en un temps minimum

#Calcul des coefficient (p,v)
def coef_p(p,v):
    return (2*f*m[v]*N[p,v]/(n[p]*(a[p]**2)))
    
#Calcul des coefficient [p,v]
def coef_c(p,v):
    return (2*f*m[v]*P[p,v]/(n[p]*(a[p]**2)))

#Calcul de A
def A(l):
    n = len(l)
    A = np.zeros((n,n))
    for i in l:
        for j in l:
            if i==j:
                s_i = 0
                for k in l:
                    if k != i:
                        s_i += coef_p(i,k)
                A[i,i] = s_i
            else:
                A[i,j] = coef_c(i,j)
    return A

def resol_tot_H(l):
    n=len(l)
    A=A(l)
    B=np.dot(A,A)
    iP,D,P = diagonalisation(B)
    def resol_K(l):
        K = np.zeros((n,3))
        for i in range(n):
            wi = np.sqrt(D[i,i])
            K[i,2] = w1                         #pulsation
            for j in range(n):
                K[i,0] += P[i,j]*H0[l[j]]       #Ai
                Cj = 0
                for k in range(n):
                    Cj += A[j,k]*L0[l[k]]       #hj'(0)
                K[i,1] += P[i,j]*Cj
                K[i,1] = K[i,1]/wi
        return K
    def resol_J(l):
        J = np.zeros((n,3))
        for i in range(n):
            wi = np.sqrt(D[i,i])
            J[i,2] = w1                         #pulsation
            for j in range(n):
                J[i,0] += P[i,j]*L0[l[j]]       #Ai
                Cj = 0
                for k in range(n):
                    Cj += A[j,k]*H0[l[k]]       #lj'(0)
                J[i,1] += P[i,j]*Cj
                J[i,1] = J[i,1]/wi
        return J
    K=resol_K(l)
    J=resol_J(l)
    def resol_H():          #Va renvoyer un tableau de fonctions
            H=[]
            for i in range(n):
                for j in range(n):
                    def s(t):
                        return iP[i,j]*(K[i,0]*np.cos(K[i,2]*t)+K[i,1]*np.sin(K[i,2]*t))
                    def Hj(t)=
                H.append(hi)
            return H
        