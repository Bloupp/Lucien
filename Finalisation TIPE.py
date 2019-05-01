import numpy as np
from numpy.polynomial import Polynomial
import numpy.linalg as alg
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
import scipy.integrate as integr

pi=np.pi
annee = 31557600

##Méthode de détermination des valeurs propres par l'algorithme de Leverrier :

#On fait tout d'abord une fonction permettant de calculer la suite des valeurs Mk que l'on range dans une liste.
#Cette fonction s'effectue en O(n**3)
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

def Dm(M):
    n = len(M)
    D = np.zeros((n,n))
    vp = valeurs_propres( M )
    for i in range(n):
        D[i,i] = vp[i]
    return D
   
''' Calcul de la base de vecteurs propres infructueux
def Qm(M):
    n = len(list)
    vp = valeurs_propres( M )
    X0 = np.zeros((1,n))
    Q = np.zeros((n,0))
    for i in range(n):
        S = M-vp[i]*np.eye(n)
        X = np.linalg.solve(S,X0)
        Q = np.concatenate((Q,X),axis=1)
    return Q
'''

def diagonalisation(M):
    D,Q=alg.eig(M)
    return alg.inv(Q),np.diag(D),Q
    

## Je fais une fonction qui va renvoyer la matrice colonne des fonctions de H et L

#La fonction qui suit se place comme un moyen d'application de l'algorithme construit jusqu'ici, celle-ci est faite dans le but de s'executer en un temps minimum

#On a tout d'abord besoin des valeurs
G = 6.67408*(10**-11)
Ms = 1.989*(10**30)
m = [3.285*(10**23), 4.867*(10**24), 5.972*(10**24), 6.39*(10**23), 1.898*(10**27), 5.863*(10**26), 8.681*(10**25), 1.024*(10**26)]
a = [57909227000, 108208475000, 149598262000, 227943824000, 778340821000, 1426666422000, 2870658186000, 4498396441000]
n_p = [ np.sqrt(G*Ms/(a[i]**3)) for i in range(8) ]
T_j = [87.95569, 224.667, 365.256363, 686.885, 4332.01, 10754, 30698 ,60216.8 ]
T_s = [24*3600*T for T in T_j]
n_p2= [ 2*pi/T for T in T_s]

e = [0.20564, 0.0068, 0.0167, 0.0934, 0.0484, 0.0538, 0.0472, 0.0086]
lpd = [77.43, 131.6, 102.937, 336.1, 14.755, 92.64, 170.92, 44.984]
lpr = [x*pi/180 for x in lpd]

H0 = [e[i]*np.sin(lpr[i]) for i in range(8)]
L0 = [e[i]*np.cos(lpr[i]) for i in range(8)]

planetes = [0,1,2,3,4,5,6,7]

#Calcul des coefficients de Fourier, utiles lors du calcul des coefficients :
def Fourier(n,i,j):
    alp = a[i]/a[j]
    def signal(t):
        return alp*np.cos(n*t)/(a[j]*((1+(alp**2)-2*alp*np.cos(t))**(3/2)))
    return integr.quad(signal,-np.pi,np.pi)[0]/np.pi

#Calcul du tableau de coefficients N(p,v)
def N():
    n = len(m)
    N=np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            Nij = Fourier(1,i,j)/16
            N[i,j],N[j,i] = Nij, Nij
    for i in range(n):
        N[i,i] = 0
    return N

NL = N()

#Calcul du tableau de coefficients P(p,v)
def P():
    n = len(m)
    P=np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            Pij = Fourier(2,i,j)/16
            P[i,j],P[j,i] = Pij, Pij
    for i in range(n):
        P[i,i] = 0
    return P
    
PL = P()

#Calcul des coefficient (p,v)
def coef_p(p,v):
    return (2*G*m[v]*NL[p,v]/(n_p[p]*(a[p]**2)))
    
#Calcul des coefficient [p,v]
def coef_c(p,v):
    return (2*G*m[v]*PL[p,v]/(n_p[p]*(a[p]**2)))

#Calcul de A avec le tableur de planètes prises en compte
def Am(l):
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
                A[i,j] = -coef_c(i,j)
    return A

def add_func(f,g):
    return lambda x: f(x) + g(x)

#Je refait la fonction principale, une fois toutes les fonctions annexes finies, celle-ci prend en entrée une liste de planètes à prendre en compte et renvoie deux listes de fonctions H(t) et L(t)
def solution_HL(l):
    n = len(planetes)
    A = Am(planetes)
    iQ,D,Q = diagonalisation(A)
    #On construit les conditions initiales :
    Y0 = np.dot(iQ,H0)
    Yp0 = np.dot(np.dot(iQ,A),L0)
    Z0 = np.dot(iQ,L0)
    Zp0 = -np.dot(np.dot(iQ,A),H0)
    #On en déduit les valeurs des constantes d'intégration :
    alpha_h = Y0
    beta_h = [Yp0[i]/D[i,i] for i in range( n )]
    alpha_l = Z0
    beta_l = [Zp0[i]/D[i,i] for i in range( n )]
    #On obtient alors les tableaux de fonctions suivants
    Y = [lambda t : alpha_h[i]*np.cos(D[i,i]*t) + beta_h[i] *np.sin(D[i,i]*t) for i in range(n)]
    Z = [lambda t : alpha_l[i]*np.cos(D[i,i]*t) + beta_l[i] *np.sin(D[i,i]*t)for i in range(n)]
    #On recombine les différentes valeurs pour obtenir H(t) = QY(t), et L(t) = QZ(t)
    H,L = [lambda x : 0 for i in range(n)],[lambda x : 0 for i in range(n)]
    for i in range(n):
        for j in range(n):
            H[i] = add_func(H[i],lambda x : Q[i,j]*Y[j](x))
            L[i] = add_func(L[i],lambda x : Q[i,j]*Z[j](x))
    return H,L

H,L = solution_HL(planetes)
T=np.linspace(-100000*annee,100000*annee,20000)     

plt.figure()
E2 = [np.sqrt(H[2](t)**2 + L[2](t)**2) for t in T]
plt.plot(T,E2, lw=2)
#plt.plot(T,E3, lw=2)
plt.show()

