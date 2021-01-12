from math import *

class Point:
    Lat, Long, Az, Beta1, Beta0, W2, A1, B1, Sigma, alphaE, Sigma1, Sm, Dsigma, Sigmaf, x, y, z = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
def geodesic(lambda1, phi1, alpha1, s, a, b):
    f = (a-b)/a
    ep = sqrt((a**2-b**2)/b**2)
    eps = 0.01
    beta1 = atan((1-f)*tan(phi1))
    beta0 = acos(cos(beta1)*sin(alpha1))
    w2 = (ep*sin(beta0))**2
    sigma1 = atan(tan(beta1)/cos(alpha1))
    alphae = asin(cos(beta0))
    Ap = 1 + w2*(4096+w2*(-768+w2*(320-175*w2)))/16384
    Bp = w2*(256+w2*(-128+w2*(74-47*w2)))/1024
    sigma = s/(b*Ap)
    sigmam = (2*sigma1+sigma)/2
    T = []
    delsig = Bp*sin(sigma)*(cos(2*sigmam)+0.25*Bp*(cos(sigma)*(2*cos(2*sigmam)**2-1)-Bp/6*cos(2*sigmam)*(-3+4*sin(sigma)**2)*(-3+4*cos(2*sigmam)**2)))
    while delsig>eps:
        T.append(delsig)
        sigma = s/(b*Ap)+delsig
        sigmam = (2*sigma1+sigma)/2
        delsig = Bp*sin(sigma)*(cos(2*sigmam)+0.25*Bp*(cos(sigma)*(2*cos(2*sigmam)**2-1)-Bp/6*cos(2*sigmam)*(-3+4*sin(sigma)**2)*(-3+4*cos(2*sigmam)**2)))
    
    tanb2 = (sin(beta1)*cos(sigma)+cos(beta1)*sin(sigma)*cos(alpha1))/sqrt(sin(alphae)**2+(sin(beta1)*sin(sigma)-cos(beta1)*cos(sigma)*cos(alpha1))**2)

    tanphi2 = tanb2/(1-f)

    tandelu = sin(sigma)*sin(alpha1)/(cos(beta1)*cos(sigma)-sin(beta1)*sin(sigma)*cos(alpha1))

    C = f/16*cos(alphae)**2*(4+f*(4-3*cos(alphae)**2))

    dellambda = atan(tandelu) - (1-C)*f*sin(alphae)*(sigma+C*sin(sigma)*(cos(2*sigmam)+C*cos(sigma)*(-1+2*cos(2*sigmam)**2)))
    lamf =lambda1+dellambda
    if lamf<0 and phi1>0:
    	lamf=lamf+pi
    elif lamf<0 and phi1<0:
    	lamf=lamf

    alpha2 = atan(sin(alphae)/(cos(beta1)*cos(sigma)*cos(alpha1)-sin(beta1)*sin(sigma)))
    
    alpha2 = alpha2 + pi
    if ( alpha2 < 0.0 ) :
    	alpha2 = alpha2 + 2*pi
    if ( alpha2 > 2*pi ) :
    	alpha2 = alpha2 - 2*pi
    return atan(tanphi2), lamf, alpha2

def geodesicpoints(a, b, phi1, lambda1, az1, s):


    f = (a-b)/a
    #1ere eccentricitee e
    e2 = sqrt((a**2-b**2)/a**2)
    #2eme eccentricitee e'
    e1 = sqrt((a**2-b**2)/b**2)

    #coordonnees du point de depart lat long
    #distance geodesique
    #azimut du pt de depart
    #nb de pts lighanstaamlu flechantillonnage, ls pts appartenant lla geodesique
    n = 400
    # pas dechanitllonnage
    s = s/(n-1)
    # array containing the Point elements
    p = []
    #intro dlpoint lwl
    P1 = Point()
    #initialaisation dyal les coordonnes o lazimut bceux du pt de depart
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
    # tant que delta sigma kber mn la precision libghina neawdu ls calculs o kula mra n updatiw lvalues d sigma, sigma_m o delta sigma
    while abs(P1.Dsigma)>0.00001:
        P1.Sigma = P1.Sigma+P1.Dsigma
        P1.Sm = (2*P1.Sigma1+P1.Sigma)/2
        P1.Dsigma = (P1.B1)*sin(P1.Sigma)*(cos(2*P1.Sm)+(P1.B1/4)*(cos(P1.Sigma)*(2*(cos(2*P1.Sm))**2-1)-(P1.B1)/6*cos(2*P1.Sm)*(-3+4*(sin(P1.Sigma))**2)*(-3+4*(cos(2*P1.Sm))**2)))

    #la valeur finale de sigma hya akhir sigma lqinaha mn loop
    P1.Sigmaf = P1.Sigma
    #addin d Point to d array
    p.append(P1)
    #redo same for the n points again calculs flcours
    for i in range(1,n):
        Pi = Point()
        Pi.Lat, Pi.Long, Pi.Az = geodesic(p[i-1].Long, p[i-1].Lat, p[i-1].Az, s, a, b)
        

        p.append(Pi)
    


    return p