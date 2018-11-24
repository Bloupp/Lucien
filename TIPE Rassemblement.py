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

#On résoud le système. Cette fonction revoie 2 tableaux de fonctions
def resol_totale(l):
    n=len(l)
    A=A(l)
    B=np.dot(A,A)
    iP,D,P = diagonalisation(B)
    def resol_red(R1,R2):
        K = np.zeros((n,3))
        for i in range(n):
            wi = np.sqrt(D[i,i])
            K[i,2] = w1                         #pulsation
            for j in range(n):
                K[i,0] += P[i,j]*R1[l[j]]       #Ai
                Cj = 0
                for k in range(n):
                    Cj += A[j,k]*R2[l[k]]
                K[i,1] += P[i,j]*Cj
                K[i,1] = K[i,1]/wi              #Bi
        return K
    K=resol_red(H0,L0)
    J=resol_red(L0,H0)
    def dered(C):          #Va renvoyer un tableau de fonctions
        def add_func(f,g):
            return lambda x: f(x) + g(x)
        H=[lambda x:0 for i in range(n)]
        for i in range(n):
            for j in range(n):
                H[i] = add_func (H[i],(lambda x: iP[i,j]*(C[j,0]*np.cos(C[j,2]*x)+C[j,1]*np.sin(C[j,2]*x))))
        return H
    return dered(K), dered(J)

'''
H,L=resol_totale(l)

T=np.linspace(0,3.153*(10**12),10)
'''
