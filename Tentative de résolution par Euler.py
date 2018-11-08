
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