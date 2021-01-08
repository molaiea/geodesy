from math import *
import numpy as np
def geodesic(lambda1, phi1, alpha1, s, a, b):
    f = (a-b)/a
    ep = sqrt((a**2-b**2)/b**2)
    eps = 0.0001
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

    # alpha2 = atan(sin(alphae)/(cos(beta1)*cos(sigma)*cos(alpha1)-sin(beta1)*sin(sigma)))
    return atan(tanphi2), lambda1+dellambda

print(geodesic(pi/6,2*pi/9,0.2*pi,922300,6378137,6356752.314245))
