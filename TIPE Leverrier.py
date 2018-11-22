import numpy as np
from numpy.polynomial import Polynomial
import numpy.linalg as alg
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D


## Valeurs
G = 6.67E-11
MS = 1.9884E30
MT = 5.972E24
PT = 1.741E11
VP = 30.4E3

#Faire une BDD contenant les valeurs utiles au calcul de (p,v) et [p,v]
#On va avoir un tableau m,n,a,N,P contenant les valeurs

#On va avoir des tableaux contenant les CI de TOUTES les planètes
H0 = []
L0 = []

## Détermination de A :

#Calcul des coefficient (p,v)
def coef_p(p,v):
    return (2*f*m[v]*N[p,v]/(n[p]*(a[p]**2)))
    
#Calcul des coefficient [p,v]
def coef_c(p,v):
    return (2*f*m[v]*P[p,v]/(n[p]*(a[p]**2)))
 
#Construction de A, l=liste des numéros de planètes prises en compte
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

#Diagonalisation de A(l), l une liste de planètes donnée : matrice diagonale
def D(l):
    n = len(l)
    D = np.zeros((n,n))
    vp = valeurs_propres( B(l) )
    for i in range(n):
        D[i,i] = vp[i]
    return D
    
#Diagonalisation de A(l), l une liste de planètes donnée : Q
def Q(l):
    n = len(l)
    B=B(l)
    vp = valeurs_propres( B )
    X0 = np.zeros((1,n))
    Q = np.zeros((n,0))
    for i in range(n):
        S = B-vp[i]*np.eye(n)
        X = np.linalg.solve(S,X0)
        Q = np.concatenate(Q,X,axis=1)
    return Q

#Résolution de K : Cette fonction va renvoyer un tableau contenant (Ai,Bi,wi) à la ieme ligne, les constantes de résolution. La connaissance de celà nous permettera de connaitre K entièrement et donc de connaitre H entièrement
def resol_K(l):
    n = len(l)
    A = A(l)
    D = D(l)
    P = Q(l)
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
    
#Résolution de H : H=P^{-1]*K  --> A optimiser
#Cette fonction renvoie un tableau où la ieme ligne est la valeur prise par les h(l(i)) à l'instant t 
def résol_H(l,t):
    n=len(l)
    K = resol_K(l)
    inv_P = alg.inv(Q(l))
    H = []
    for i in range(n):
        hi=0
        for j in range(n):
            hi += inv_P[i,j]*(K[i,0]*np.cos(K[i,2]*t)+K[i,1]*np.sin(K[i,2]*t))
        H.append(hi)
    return H
    
#On fait pareil pour L, la seule différence va être les conditions initiales vue que les équations différentielles sont les mêmes
def resol_J(l):
    n = len(l)
    A = A(l)
    D = D(l)
    P = Q(l)
    J = np.zeros((n,3))
    for i in range(n):
        wi = np.sqrt(D[i,i])
        J[i,2] = w1                         #pulsation
        for j in range(n):
            J[i,0] += P[i,j]*L0[l[j]]       #Ai
            Cj = 0
            for k in range(n):
                Cj += A[j,k]*H0[l[k]]       #hj'(0)
            J[i,1] += P[i,j]*Cj
            J[i,1] = J[i,1]/wi
    return K
    
#Résolution de L : L=P^{-1]*J  --> A optimiser
#Cette fonction renvoie un tableau où la ieme ligne est la valeur prise par les L(l(i)) à l'instant t 
def résol_L(l,t):
    n=len(l)
    J = resol_J(l)
    inv_P = alg.inv(Q(l))
    L = []
    for i in range(n):
        hi=0
        for j in range(n):
            li += inv_P[i,j]*(J[i,0]*np.cos(J[i,2]*t)+J[i,1]*np.sin(J[i,2]*t))
        L.append(li)
    return L
