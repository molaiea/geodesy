from math import *
from . import directp

def geodesic(lambda1, phi1, alpha1, s, a, b):
    f = (a-b)/a
    ep = sqrt((a**2-b**2)/b**2)
    eps = 0.01
    #Calcul de la latitude réduite beta1 
    beta1 = atan((1-f)*tan(phi1))
    #Calcul de la latitude réduite beta0
    beta0 = acos(cos(beta1)*sin(alpha1))
     #Calcul de la constante w2 
    w2 = (ep*sin(beta0))**2
    #Calcul de la distance angulaire sigma1 sur la sphère auxiliaire 
    sigma1 = atan(tan(beta1)/cos(alpha1))
    #Calcul de l'azimut de la géodésique à l'équateur alphae
    alphae = asin(cos(beta0))
    #Calcul des constantes de Vincenty Ap et Bp
    Ap = 1 + w2*(4096+w2*(-768+w2*(320-175*w2)))/16384
    Bp = w2*(256+w2*(-128+w2*(74-47*w2)))/1024
    #Calcul de la distance angulaire sigma sur le grand cercle de la sphère auxiliaire 
    sigma = s/(b*Ap)
    sigmam = (2*sigma1+sigma)/2
    delsig = Bp*sin(sigma)*(cos(2*sigmam)+0.25*Bp*(cos(sigma)*(2*cos(2*sigmam)**2-1)-Bp/6*cos(2*sigmam)*(-3+4*sin(sigma)**2)*(-3+4*cos(2*sigmam)**2)))
    while delsig>eps:
        sigma = s/(b*Ap)+delsig
        sigmam = (2*sigma1+sigma)/2
        delsig = Bp*sin(sigma)*(cos(2*sigmam)+0.25*Bp*(cos(sigma)*(2*cos(2*sigmam)**2-1)-Bp/6*cos(2*sigmam)*(-3+4*sin(sigma)**2)*(-3+4*cos(2*sigmam)**2)))
    #Calcul de tanb2
    tanb2 = (sin(beta1)*cos(sigma)+cos(beta1)*sin(sigma)*cos(alpha1))/sqrt(sin(alphae)**2+(sin(beta1)*sin(sigma)-cos(beta1)*cos(sigma)*cos(alpha1))**2)
    #Calcul de phi2
    tanphi2 = tanb2/(1-f)
    #Calcul de la tangente de la différence de longitude sur la sphère auxiliaire tandelu
    tandelu = sin(sigma)*sin(alpha1)/(cos(beta1)*cos(sigma)-sin(beta1)*sin(sigma)*cos(alpha1))
    #Calcul de la constante C de Vincenty 
    C = f/16*cos(alphae)**2*(4+f*(4-3*cos(alphae)**2))
    #Calcul de la différence de longitude dellambda
    dellambda = atan(tandelu) - (1-C)*f*sin(alphae)*(sigma+C*sin(sigma)*(cos(2*sigmam)+C*cos(sigma)*(-1+2*cos(2*sigmam)**2)))
    lamf =lambda1+dellambda
    if lamf<0 and phi1>0:
        lamf=lamf+pi
    elif lamf<0 and phi1<0:
        lamf=lamf
    #Calcul de l'azimut alpha2 et de l'azimut inverse 
    alpha2 = atan(sin(alphae)/(cos(beta1)*cos(sigma)*cos(alpha1)-sin(beta1)*sin(sigma)))
    az = alpha2
    if az<0 and lamf>0:
        az = az+pi
    elif az<0 and lamf<0:
        az = az-pi
    elif az>0 and lamf>0:
        az = az+pi
    else:
        az=az-pi
    return atan(tanphi2), lamf, az

class Point:
    Lat, Long, Az, Beta1, Beta0, W2, A1, B1, Sigma, alphaE, Sigma1, Sm, Dsigma, Sigmaf, x, y, z = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0


def geodesicpoints(a, b, phi1, lambda1, az1, s):


    f = (a-b)/a
    #1ere eccentricitee e
    e2 = sqrt((a**2-b**2)/a**2)
    #2eme excentricité e'
    e1 = sqrt((a**2-b**2)/b**2)
    #coordonnées du point de départ lat long
    #distance géodesique
    #azimut du pt de depart
    #nb de pts pour l'echantillonnage, ls pts appartenant à la geodesique
    n = 400
    # pas d'echantillonnage
    s = s/(n-1)
    # array contenant les elements du point
    p = []
    #intro du premier point
    P1 = Point()
    #initialaisation  des coordonnées o lazimut par ceux du du point départ
    P1.Lat = phi1
    P1.Long = lambda1
    P1.Az = az1
    #calculs necessaires pr calculer le sigma (cours)
    P1.Beta1 = atan((1-f)*tan(P1.Lat))
    P1.Beta0 = acos(cos(P1.Beta1)*sin(az1))
    P1.W2 = (e1**2)*((sin(P1.Beta0))**2)
    P1.Sigma1 = atan(tan(P1.Beta1)/cos(az1))
    P1.alphaE = asin(cos(P1.Beta0))
    P1.A1 = 1+(P1.W2/16384)*(4096+P1.W2*(-768+P1.W2*(320-175*P1.W2)))
    P1.B1 = (P1.W2/1024)*(256+P1.W2*(-128+P1.W2*(74-47*P1.W2)))
    P1.Sigma = s/(b*P1.A1)
    P1.Sm = (2*P1.Sigma1+P1.Sigma)/2
    P1.Dsigma = (P1.B1)*sin(P1.Sigma)*(cos(2*P1.Sm)+(P1.B1/4)*(cos(P1.Sigma)*(2*(cos(2*P1.Sm))**2-1)-(P1.B1)/6*cos(2*P1.Sm)*(-3+4*(sin(P1.Sigma))**2)*(-3+4*(cos(2*P1.Sm))**2)))
    # on arrete les itérations car la variation sigma devient négligeable
    while abs(P1.Dsigma)>0.001:
        P1.Sigma = P1.Sigma+P1.Dsigma
        P1.Sm = (2*P1.Sigma1+P1.Sigma)/2
        P1.Dsigma = (P1.B1)*sin(P1.Sigma)*(cos(2*P1.Sm)+(P1.B1/4)*(cos(P1.Sigma)*(2*(cos(2*P1.Sm))**2-1)-(P1.B1)/6*cos(2*P1.Sm)*(-3+4*(sin(P1.Sigma))**2)*(-3+4*(cos(2*P1.Sm))**2)))

    #la valeur finale de sigma celle qui nous donne la boucle 
    P1.Sigmaf = P1.Sigma
    #addition du Point  dans le tableau
    p.append(P1)
    #De meme pour les n points restants
    if P1.Lat>0:
        for i in range(1,n):
            lat, long, az = geodesic(p[i-1].Long, p[i-1].Lat, p[i-1].Az, s, a, b)
            Pi = Point()
            Pi.Long, Pi.Lat, Pi.Az = long, lat, az
            p.append(Pi)
    else:
        for i in range(1,n):
        # lat, long, az = geodesic(p[i-1].Long, p[i-1].Lat, p[i-1].Az, s, a, b)
            Pi = Point()
            # Pi.Lat, Pi.Long, Pi.Az = lat, long, az
            Pi = Point()
            Pi.Beta1 = atan((sin(p[i-1].Beta1)*cos(p[i-1].Sigmaf)+cos(p[i-1].Beta1)*sin(p[i-1].Sigmaf)*cos(p[i-1].Az))/sqrt((sin(p[i-1].alphaE))**2+(sin(p[i-1].Beta1)*sin(p[i-1].Sigmaf)-cos(p[i-1].Beta1)*cos(p[i-1].Sigmaf)*cos(p[i-1].Az))**2))
            Pi.Lat = atan(tan(Pi.Beta1)/(1-f))

            du = atan((sin(p[i-1].Sigmaf)*sin(p[i-1].Az))/(cos(p[i-1].Beta1)*cos(p[i-1].Sigmaf)-sin(p[i-1].Beta1)*sin(p[i-1].Sigmaf)*cos(p[i-1].Az)))

            c = f*((cos(p[i-1].alphaE))**2)/16*(4+f*(4-3*(cos(p[i-1].alphaE))**2))

            dlambda = du - (1-c)*f*sin(p[i-1].alphaE)*(p[i-1].Sigmaf+c*sin(p[i-1].Sigmaf)*(cos(2*p[i-1].Sm)+c*cos(p[i-1].Sigmaf)*(-1+2*(cos(2*p[i-1].Sm))**2)))

            Pi.Long = dlambda+p[i-1].Long
            Pi.Az = atan(sin(p[i-1].alphaE)/(cos(p[i-1].Beta1)*cos(p[i-1].Sigmaf)*cos(p[i-1].Az)-sin(p[i-1].Beta1)*sin(p[i-1].Sigmaf)))
            Pi.Beta0 = acos(cos(Pi.Beta1)*sin(Pi.Az))
            Pi.W2 = e1**2*(sin(Pi.Beta1)**2)
            Pi.Sigma1 = atan(tan(Pi.Beta1)/cos(Pi.Az))
            Pi.alphaE = asin(cos(Pi.Beta0))
            Pi.A1 = 1+(Pi.W2/16384)*(4096+Pi.W2*(-768+Pi.W2*(320-175*Pi.W2)))
            Pi.B1 = (Pi.W2/1024)*(256+Pi.W2*(-128+Pi.W2*(74-47*Pi.W2)))
            #sigma la distance angulaire entre 2 points consecutives est la meme s/(b*Pi.A1) s etant le pas
            Pi.Sigma = s/(b*Pi.A1)
            Pi.Sm = (2*Pi.Sigma1+Pi.Sigma)/2
            Pi.Dsigma = (Pi.B1)*sin(Pi.Sigma)*(cos(2*Pi.Sm)+(Pi.B1/4)*(cos(Pi.Sigma)*(2*(cos(2*Pi.Sm))**2-1)-(Pi.B1)/6*cos(2*Pi.Sm)*(-3+4*(sin(Pi.Sigma))**2)*(-3+4*(cos(2*Pi.Sm))**2)))
            while abs(Pi.Dsigma)>0.00001:
                Pi.Sigma = Pi.Sigma+Pi.Dsigma
                Pi.Sm = (2*Pi.Sigma1+Pi.Sigma)/2
                Pi.Dsigma = (Pi.B1)*sin(Pi.Sigma)*(cos(2*Pi.Sm)+(Pi.B1/4)*(cos(Pi.Sigma)*(2*(cos(2*Pi.Sm))**2-1)-(Pi.B1)/6*cos(2*Pi.Sm)*(-3+4*(sin(Pi.Sigma))**2)*(-3+4*(cos(2*Pi.Sm))**2)))

            Pi.Sigmaf = Pi.Sigma
            if(Pi.Long)<0:
                Pi.Long=Pi.Long+pi
            p.append(Pi)
    


    return p