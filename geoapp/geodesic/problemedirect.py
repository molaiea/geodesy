from math import *
from . import direct, directp

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
    delsig = Bp*sin(sigma)*(cos(2*sigmam)+0.25*Bp*(cos(sigma)*(2*cos(2*sigmam)**2-1)-Bp/6*cos(2*sigmam)*(-3+4*sin(sigma)**2)*(-3+4*cos(2*sigmam)**2)))
    while delsig>eps:
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
    while abs(P1.Dsigma)>0.001:
        P1.Sigma = P1.Sigma+P1.Dsigma
        P1.Sm = (2*P1.Sigma1+P1.Sigma)/2
        P1.Dsigma = (P1.B1)*sin(P1.Sigma)*(cos(2*P1.Sm)+(P1.B1/4)*(cos(P1.Sigma)*(2*(cos(2*P1.Sm))**2-1)-(P1.B1)/6*cos(2*P1.Sm)*(-3+4*(sin(P1.Sigma))**2)*(-3+4*(cos(2*P1.Sm))**2)))

    #la valeur finale de sigma hya akhir sigma lqinaha mn loop
    P1.Sigmaf = P1.Sigma
    #addin d Point to d array
    p.append(P1)
    #redo same for the n points again calculs flcours
    if P1.Lat>0:
        for i in range(1,n):
            lat, long, az = geodesic(p[i-1].Long, p[i-1].Lat, p[i-1].Az, s, a, b)
            Pi = Point()
            # Pi.Long, Pi.Lat, Pi.Az = directp.direct(a,b,p[i-1].Lat, p[i-1].Long, p[i-1].Az, s)
            # Pi.Long = Pi.Long*pi/180
            # Pi.Lat = Pi.Lat*pi/180
            # Pi.Az = Pi.Az*pi/180
            # Pi = Point()
            Pi.Long, Pi.Lat, Pi.Az = long, lat, az
            # Pi.Beta1 = atan((sin(p[i-1].Beta1)*cos(p[i-1].Sigmaf)+cos(p[i-1].Beta1)*sin(p[i-1].Sigmaf)*cos(p[i-1].Az))/sqrt((sin(p[i-1].alphaE))**2+(sin(p[i-1].Beta1)*sin(p[i-1].Sigmaf)-cos(p[i-1].Beta1)*cos(p[i-1].Sigmaf)*cos(p[i-1].Az))**2))
            # Pi.Lat = atan(tan(Pi.Beta1)/(1-f))

            # du = atan((sin(p[i-1].Sigmaf)*sin(p[i-1].Az))/(cos(p[i-1].Beta1)*cos(p[i-1].Sigmaf)-sin(p[i-1].Beta1)*sin(p[i-1].Sigmaf)*cos(p[i-1].Az)))

            # c = f*((cos(p[i-1].alphaE))**2)/16*(4+f*(4-3*(cos(p[i-1].alphaE))**2))

            # dlambda = du - (1-c)*f*sin(p[i-1].alphaE)*(p[i-1].Sigmaf+c*sin(p[i-1].Sigmaf)*(cos(2*p[i-1].Sm)+c*cos(p[i-1].Sigmaf)*(-1+2*(cos(2*p[i-1].Sm))**2)))

            # Pi.Long = dlambda+p[i-1].Long
            # Pi.Az = atan(sin(p[i-1].alphaE)/(cos(p[i-1].Beta1)*cos(p[i-1].Sigmaf)*cos(p[i-1].Az)-sin(p[i-1].Beta1)*sin(p[i-1].Sigmaf)))


            # Pi.Beta0 = acos(cos(Pi.Beta1)*sin(Pi.Az))
            # Pi.W2 = e1**2*(sin(Pi.Beta1)**2)
            # Pi.Sigma1 = atan(tan(Pi.Beta1)/cos(Pi.Az))
            # Pi.alphaE = asin(cos(Pi.Beta0))
            # Pi.A1 = 1+(Pi.W2/16384)*(4096+Pi.W2*(-768+Pi.W2*(320-175*Pi.W2)))
            # Pi.B1 = (Pi.W2/1024)*(256+Pi.W2*(-128+Pi.W2*(74-47*Pi.W2)))
            # #sigma la distance angulaire between 2 consecutive points hya dima s/(b*Pi.A1) s etant le pas, ls pts dyulna aykunu equidistants
            # Pi.Sigma = s/(b*Pi.A1)
            # Pi.Sm = (2*Pi.Sigma1+Pi.Sigma)/2
            # Pi.Dsigma = (Pi.B1)*sin(Pi.Sigma)*(cos(2*Pi.Sm)+(Pi.B1/4)*(cos(Pi.Sigma)*(2*(cos(2*Pi.Sm))**2-1)-(Pi.B1)/6*cos(2*Pi.Sm)*(-3+4*(sin(Pi.Sigma))**2)*(-3+4*(cos(2*Pi.Sm))**2)))
            # while abs(Pi.Dsigma)>0.00001:
            #     Pi.Sigma = Pi.Sigma+Pi.Dsigma
            #     Pi.Sm = (2*Pi.Sigma1+Pi.Sigma)/2
            #     Pi.Dsigma = (Pi.B1)*sin(Pi.Sigma)*(cos(2*Pi.Sm)+(Pi.B1/4)*(cos(Pi.Sigma)*(2*(cos(2*Pi.Sm))**2-1)-(Pi.B1)/6*cos(2*Pi.Sm)*(-3+4*(sin(Pi.Sigma))**2)*(-3+4*(cos(2*Pi.Sm))**2)))

            # Pi.Sigmaf = Pi.Sigma
            # if(Pi.Long)<0:
            #     Pi.Long=Pi.Long+pi
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
            #sigma la distance angulaire between 2 consecutive points hya dima s/(b*Pi.A1) s etant le pas, ls pts dyulna aykunu equidistants
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