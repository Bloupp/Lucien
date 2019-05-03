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
   
def Qm(M):      #Calcul de la base de vecteurs propres infructueux
    n = len(list)
    vp = valeurs_propres( M )
    X0 = np.zeros((1,n))
    Q = np.zeros((n,0))
    for i in range(n):
        S = M-vp[i]*np.eye(n)
        X = np.linalg.solve(S,X0)
        Q = np.concatenate((Q,X),axis=1)
    return Q

def diagonalisation(M):
    D,Q=alg.eig(M)
    return alg.inv(Q),np.diag(D),Q
    


## Quelques valeurs utiles pour le système solaire
G = 6.67408*(10**-11)
Ms = 1.989*(10**30)
m = [3.285*(10**23), 4.867*(10**24), 5.972*(10**24), 6.39*(10**23), 1.898*(10**27), 5.863*(10**26), 8.681*(10**25), 1.024*(10**26)]
a = [57909227000, 108208475000, 149598262000, 227943824000, 778340821000, 1426666422000, 2870658186000, 4498396441000]
n_p = [ np.sqrt(G*Ms/(a[i]**3)) for i in range(8) ]
''' Je vérifie ici les valeurs de n_p, cela convient bien '''
T_j = [87.95569, 224.667, 365.256363, 686.885, 4332.01, 10754, 30698 ,60216.8 ]
T_s = [24*3600*T for T in T_j]
n_p2= [ 2*pi/T for T in T_s]

e = [0.20564, 0.0068, 0.0167, 0.0934, 0.0484, 0.0538, 0.0472, 0.0086]
lpd = [77.43, 131.6, 102.937, 336.1, 14.755, 92.64, 170.92, 44.984]
lpr = [x*pi/180 for x in lpd]

H0 = [e[i]*np.sin(lpr[i]) for i in range(8)]
L0 = [e[i]*np.cos(lpr[i]) for i in range(8)]



## Fonction de résolution du problème posé

#La fonction qui suit se place comme un moyen d'application de l'algorithme construit jusqu'ici, celle-ci est faite dans le but de s'executer en un temps minimum

#Calcul des coefficients de Fourier, utiles lors du calcul des coefficients :
def integrale(f,n,a,b):
    dt = (b-a)/n
    s=0
    for i in range(n):
        s += (f((i+1)*dt)+f(i*dt))*dt/2
    return s

def Fourier(n,i,j):
    alp = a[i]/a[j]
    def s(t):
        return np.cos(n*t)/((1+alp**2-2*alp*np.cos(t))**(3/2))
    cn = 2*(integrale(s,10000,0,0.1)+integrale(s,1000,0.1,pi))/pi
    return alp*cn/a[j]

#Calcul du tableau de coefficients N(p,v)
def N():
    n = len(m)
    N=np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            Nij = Fourier(1,i,j)/8
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
            Pij = Fourier(2,i,j)/8
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

def mapp(f,T):
    return [f(t) for t in T]

#Je réalise alors la fonction principale de résolution pour déterminer H et L, celle ci renvoie 2 tableaux de fonctions
def solution_HL(l):
    n = len(l)
    A = Am(l)
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
    def h(i):
        return (lambda t : sum([Q[i,j]*(alpha_h[j]*np.cos(D[j,j]*t) + beta_h[j] *np.sin(D[j,j]*t)) for j in range(n)] ))
    def l(i):
        return (lambda t : sum([Q[i,j]*(alpha_l[j]*np.cos(D[j,j]*t) + beta_l[j] *np.sin(D[j,j]*t)) for j in range(n)] ))
    Hsol = mapp(h,range(n))
    Lsol = mapp(l,range(n))
    return Hsol,Lsol

#La méthode est la même pour déterminer P et Q, la seule différence est dans les condition initiales
def solution_PQ(l):
    n = len(l)
    A = Am(l)
    iQ,D,Q = diagonalisation(A)
    #On construit les conditions initiales :
    Y0 = np.dot(iQ,P0)
    Yp0 = np.dot(iQ,np.dot(A,Q0))
    Z0 = np.dot(iQ,Q0)
    Zp0 = -np.dot(iQ,np.dot(A,P0))
    #On en déduit les valeurs des constantes d'intégration :
    alpha_h = Y0
    beta_h = [Yp0[i]/D[i,i] for i in range( n )]
    alpha_l = Z0
    beta_l = [Zp0[i]/D[i,i] for i in range( n )]
    #On obtient alors les tableaux de fonctions suivants
    def h(i):
        return (lambda t : sum([Q[i,j]*(alpha_h[j]*np.cos(D[j,j]*t) + beta_h[j] *np.sin(D[j,j]*t)) for j in range(n)] ))
    def l(i):
        return (lambda t : sum([Q[i,j]*(alpha_l[j]*np.cos(D[j,j]*t) + beta_l[j] *np.sin(D[j,j]*t)) for j in range(n)] ))
    Psol = mapp(h,range(n))
    Qsol = mapp(l,range(n))
    return Psol,Qsol


## Calcul des tableaux de fonction résultats

planetes = [0,1,2,3,4,5,6,7]

H,L = solution_HL(planetes)
P,Q = solution_PQ(planetes)


## Affichage des réultats : tracé courbes et ellipses
T=np.linspace(-100000*annee,100000*annee,20000)     

#Affichage de toutes les excentricités sur la même courbe
plt.figure()
for i in range(len(planetes)):
    E = [np.sqrt(H[i](t)**2 + L[i](t)**2) for t in T]
    plt.plot(T,E, lw=2)
plt.grid()
plt.show()


## Résultats de Le Verrier :
E_terre=[0.046,0.0473,0.0452,0.0398,0.0316,0.0218,0.0131,0.0109,0.0151,0.0188,0.0187,0.0168,0.0115,0.0047,0.0059,0.0124,0.0173,0.0199,0.0211,0.0188,0.0176,0.0189]
X=np.linspace(-110000,100000,22)

plt.figure()
plt.plot(X,E_terre)
plt.show()
