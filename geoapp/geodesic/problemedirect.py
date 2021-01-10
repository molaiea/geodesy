from math import *

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