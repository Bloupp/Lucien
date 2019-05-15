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

G = 6.67408*(10**-11)   #Constante de gravitation universelle
Ms = 1.989*(10**30)     #Masse du soleil
nom = ['Mercure', 'Vénus', 'Terre', 'Mars', 'Jupiter', 'Saturne', 'Uranus', 'Neptune']
#Masses
m = [3.285*(10**23), 4.867*(10**24), 5.972*(10**24), 6.39*(10**23), 1.898*(10**27), 5.863*(10**26), 8.681*(10**25), 1.024*(10**26)]
#Demi-grands axes
a = [57909227000, 108208475000, 149598262000, 227943824000, 778340821000, 1426666422000, 2870658186000, 4498396441000]
#Moyen mouvements
n_p = [ np.sqrt(G*Ms/(a[i]**3)) for i in range(8) ]
''' Je vérifie ici les valeurs de n_p, cela convient bien'''
#Périodes de révolution en jours
T_j = [87.95569, 224.667, 365.256363, 686.885, 4332.01, 10754, 30698 ,60216.8 ]
#Périodes de révolution en secondes
T_s = [24*3600*T for T in T_j]
#Moyen mouvements
n_p2= [ 2*pi/T for T in T_s]

#Excentricités
e = [0.20564, 0.0068, 0.0167, 0.0934, 0.0484, 0.0538, 0.0472, 0.0086]
#Longitudes du périhélie en degrés
lpd = [77.43, 131.6, 102.937, 336.1, 14.755, 92.64, 170.92, 44.984]
#Longitudes du périhélie en radians
lpr = [x*pi/180 for x in lpd]

#Valeurs initiales de H et L
H0 = [e[i]*np.sin(lpr[i]) for i in range(8)]
L0 = [e[i]*np.cos(lpr[i]) for i in range(8)]

#Plan de l'obite en degrés (Noeud ascendant)
thd = [48.33, 76.7, 174.873, 49.6, 100.5, 113.7, 74.02, 131.784]
#Inclianaisons de l'orbite en degrés
incld = [7.005, 3.3947, 0, 1.85, 1.304, 2.486, 0.773, 1.77]
#Plan de l'obite en radians (Noeud ascendant)
th = [pi*thd[i]/180 for i in range(8)]
#Inclianaisons de l'orbite en radians
incl = [pi*incld[i] for i in range(8)]

#Valeurs initiales de P et Q
P0 = [np.tan(incl[i])*np.sin(th[i]) for i in range(8)]
Q0 = [np.tan(incl[i])*np.cos(th[i]) for i in range(8)]


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
    for i in range(n):
        for j in range(n):
            if i==j:
                s_i = 0
                for k in l:
                    if k != i:
                        s_i += coef_p(l[i],k)
                A[i,i] = s_i
            else:
                A[i,j] = -coef_c(l[i],l[j])
    return A

def mapp(f,T):
    return [f(t) for t in T]

#Je réalise alors la fonction principale de résolution pour déterminer H et L, celle ci renvoie 2 tableaux de fonctions
#Ici CI1 et CI2 sont les conditions initiales
def solution(l,CI1,CI2):
    n = len(l)
    #On construit le tableau de constantes convenant puis on découpe les conditions initiales a la bonne taille
    A = Am(l)
    iQ,D,Q = diagonalisation(A)
    C1 , C2 = [CI1[i] for i in l], [CI2[i] for i in l]
    #On construit les conditions initiales :
    Y0 = np.dot(iQ,C1)
    Yp0 = np.dot(np.dot(iQ,A),C2)
    Z0 = np.dot(iQ,C2)
    Zp0 = -np.dot(np.dot(iQ,A),C1)
    #On en déduit les valeurs des constantes d'intégration :
    alpha_h = Y0
    beta_h = [Yp0[i]/D[i,i] for i in range( n )]
    alpha_l = Z0
    beta_l = [Zp0[i]/D[i,i] for i in range( n )]
    #On obtient alors les tableaux de fonctions suivants
    Y = [lambda t : alpha_h[i]*np.cos(D[i,i]*t) + beta_h[i] *np.sin(D[i,i]*t) for i in range(n)]
    Z = [lambda t : alpha_l[i]*np.cos(D[i,i]*t) + beta_l[i] *np.sin(D[i,i]*t)for i in range(n)]
    #On recombine les différentes valeurs pour obtenir H(t) = QY(t), et L(t) = QZ(t)
    def sol1(i):
        return (lambda t : sum([Q[i,j]*(alpha_h[j]*np.cos(D[j,j]*t) + beta_h[j] *np.sin(D[j,j]*t)) for j in range(n)] ))
    def sol2(i):
        return (lambda t : sum([Q[i,j]*(alpha_l[j]*np.cos(D[j,j]*t) + beta_l[j] *np.sin(D[j,j]*t)) for j in range(n)] ))
    Sol1 = mapp(sol1,range(n))
    Sol2 = mapp(sol2,range(n))
    return Sol1,Sol2

#La méthode est la même pour déterminer P et Q, la seule différence est dans les condition initiales

## Calcul des tableaux de fonction résultats

planetes = [0,1,2,3,4,5,6,7]

H,L = solution(planetes,H0,L0)
P,Q = solution(planetes,P0,Q0)

def excentricite(i,t,R):
    return np.sqrt(R[0][i](t)**2 + R[1][i](t)**2)

#Les fonctions suivantes donnent des résultats étranges... Il faut ici gérer les cas ou tan passe par pi/2 à la main
def longitudeperihelie(i,t,R):
    v = R[0][i](t)/R[1][i](t)
    return np.arctan(v)

def theta(i,t,R):
    v = R[2][i](t)/R[3][i](t)
    return np.arctan(v)

def inclinaison(i,t,R):
    v = np.sqrt(R[2][i](t)**2 + R[3][i](t)**2)
    return np.arctan(v)


## Affichage des réultats : tracé courbes et ellipses pour la configuration présente

duree = 100000
T = np.linspace(-duree*annee,duree*annee,20000)
T_aff = np.linspace(-duree,duree,20000)

#Affichage de toutes les excentricités sur la même courbe
plt.figure()
for i in range(1,len(planetes)):
    E = [excentricite(i,t,[H,L,P,Q]) for t in T]
    plt.plot(T_aff,E, lw=2)
plt.ylabel("Excentricités")
plt.xlabel("Temps - à 0 à l'instant présent (en années)")
plt.title("Variation des excentricités au cours du temps")
plt.legend(nom[1:])
plt.grid()
plt.show()

#Affichage des Noeuds ascendants
plt.figure()
for i in range(1,len(planetes)):
    N = [theta(i,t,[H,L,P,Q])*180/pi for t in T]
    plt.plot(T_aff,N, lw=2)
plt.grid()
plt.show()

#Affichage des inclinaisons
plt.figure()
for i in range(1,len(planetes)):
    I = [inclinaison(i,t,[H,L,P,Q])*180/pi for t in T]
    plt.plot(T_aff,I, lw=2)
plt.ylabel("Inclinaisons")
plt.xlabel("Temps - à 0 à l'instant présent (en années)")
plt.title("Variation des inclinaisons au cours du temps")
plt.legend(nom[1:])
plt.grid()
plt.show()

#Mettre les fonctions de tracé ellipses


## Analyse des résultats

#Pour simplifier les fonctions, on va toujours prendre en entrée le tableau des fonctions [H,L,P,Q] correspondantes où H,L,P,Q sont les tableaux de fonctions résultats pour les conditions initiales données

#En faisant uniquement varier les angles, excentricites et planetes prises en compte
def TabSol(l,H01,L01,P01,Q01):
    H1,L1 = solution(l,H01,L01)
    P1,Q1 = solution(l,P01,Q01)
    return [H1,L1,P1,Q1]

Sol1,Sol2 = TabSol(planetes,H0,L0,P0,Q0), TabSol([0,1,2,3,6,7],H0,L0,P0,Q0)

#On test si l'écart entre les excentricités de 2 solutions dépasse epsilon : i est la planete à étudier, (R1,R2) les 2 configurations considérées
def testVarExc(i,R1,R2,max,pas,epsilon):
    k = 0
    e1,e2 = excentricite(i,0,R1),excentricite(i,0,R2)
    while k>max and abs(e1-e2)<epsilon:
        e1,e2 = excentricite(i,k*pas*annee,R1),excentricite(i,k*pas*annee,R2)
        k += 1
    if k == max:
        return False
    else:
        return k

#On calcul ici l'écart moyen des excentricités de 2 solutions entre t1 et t2
def EcartMoyExc(i,R1,R2,t1,t2):
    n = abs(t2-t1)//1000
    s=0
    for j in range(n):
        s += abs(excentricite(i,j*annee*1000+t1,R1)-excentricite(i,j*annee*1000+t1,R2))
    return s/n


## Résultats de Le Verrier :
E_terre=[0.046,0.0473,0.0452,0.0398,0.0316,0.0218,0.0131,0.0109,0.0151,0.0188,0.0187,0.0168,0.0115,0.0047,0.0059,0.0124,0.0173,0.0199,0.0211,0.0188,0.0176,0.0189]
X=np.linspace(-110000,100000,22)

plt.figure()
plt.plot(X,E_terre)
plt.show()
