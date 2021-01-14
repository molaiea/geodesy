from math import *

def direct(latitude1, longitude1, alpha1, s, a, b):
    f = (a-b)/a		
    ep = sqrt((a**2-b**2)/b**2)
    #Calcul de la latitude réduite beta1 
    TanBeta1 = (1-f) * tan(latitude1)
    Beta1 = atan( TanBeta1 )
    #Calcul de la latitude réduite beta0
    CosBeta0 = cos(Beta1)*sin(alpha1)
    Beta0 = acos(CosBeta0)
    #Calcul de la constante w2 
    w2 = (sin(Beta0)*ep)**2
    #Calcul de la distance angulaire sigma1 sur la sphère auxiliaire 
    sigma1 = atan2( TanBeta1, cos(alpha1) )
    #Calcul de l'azimut de la géodésique à l'équateur alphae
    Sinalphae = CosBeta0
    alphae = asin(Sinalphae)
    cosalphae_sq = 1.0 - Sinalphae**2
    #Calcul des constantes de Vincenty A' et B'
    A = 1.0 + (w2 / 16384) * (4096 + w2 * (-768 + w2 * \
            (320 - 175 * w2) ) )
    B = (w2 / 1024) * (256 + w2 * (-128 + w2 * (74 - 47 * w2) ) )
    #Calcul de la distance angulaire sigma sur le grand cercle de la sphère auxiliaire 
    sigma = (s / (b * A))
    sigma_m = (2 * sigma1 + sigma)/2
    count = 0
    while ( abs(sigma_m) > 0.01 and count<200) :
            sigma_m = (2 * sigma1 + sigma)

            delta_sigma = B * sin(sigma) * ( cos(sigma_m) 
                    + (B/4) * (cos(sigma) * 
                    (-1 + 2 * pow( cos(sigma_m), 2 ) -  
                    (B/6) * cos(sigma_m) * \
                    (-3 + 4 * pow(sin(sigma), 2 )) *  \
                    (-3 + 4 * pow( cos (sigma_m), 2 ))))) 
            
            sigma = (s / (b * A)) + delta_sigma
            count = count+1
    #Calcul de la latitude 2
    latitude2 = atan2 ( (sin(Beta1) * cos(sigma) + cos(Beta1) * sin(sigma) * cos(alpha1) ), \
            ((1-f) * sqrt( pow(Sinalphae, 2) +  \
            pow(sin(Beta1) * sin(sigma) - cos(Beta1) * cos(sigma) * cos(alpha1), 2))))
    #Calcul de la différence de longitude sur la sphère auxiliaire deltau
    deltau = atan2( (sin(sigma) * sin(alpha1 )), (cos(Beta1) * cos(sigma) -  \
            sin(Beta1) *  sin(sigma) * cos(alpha1)))
    #Calcul de la constante C de Vincenty 
    C = (f/16) * cosalphae_sq * (4 + f * (4 - 3 * cosalphae_sq ))
    #Calcul de la différence de longitude dellambda
    deltalambda = deltau - (1-C) * f * Sinalphae *  \
            (sigma + C * sin(sigma) * (cos(sigma_m) + \
            C * cos(sigma) * (-1 + 2 * pow(cos(sigma_m),2) )))

    longitude2 = longitude1 + deltalambda
    #Calcul de l'azimut alpha2 et de l'azimut inverse 
    alpha21 = atan2 ( Sinalphae, (-sin(Beta1) * sin(sigma) +  \
            cos(Beta1) * cos(sigma) * cos(alpha1)))

    alpha21 = alpha21 + 2*pi / 2.0
    if ( alpha21 < 0.0 ) :
            alpha21 = alpha21 + 2*pi
    if ( alpha21 > 2*pi ) :
            alpha21 = alpha21 - 2*pi

    latitude2       = latitude2*180/pi
    longitude2    = longitude2    *180/pi
    alpha21    = alpha21*180/pi
    if longitude2>180:
            longitude2 = longitude2-360
    return round(longitude2), round(latitude2), round(alpha21)
