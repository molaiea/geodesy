from math import *
import numpy as np

def direct(a, b, latitude1, longitude1, alpha1To2, s ) :


        f = (a-b)/a		

        if ( alpha1To2 < 0.0 ) : 
                alpha1To2 = alpha1To2 + 2*pi
        if ( alpha1To2 > 2*pi ) : 
                alpha1To2 = alpha1To2 - 2*pi

        TanU1 = (1-f) * tan(latitude1)
        U1 = atan( TanU1 )
        sigma1 = atan2( TanU1, cos(alpha1To2) )
        Sinalpha = cos(U1) * sin(alpha1To2)
        cosalpha_sq = 1.0 - Sinalpha * Sinalpha

        u2 = cosalpha_sq * (a * a - b * b ) / (b * b)
        A = 1.0 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * \
                (320 - 175 * u2) ) )
        B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2) ) )

        # Starting with the approximation
        sigma = (s / (b * A))

        last_sigma = 2.0 * sigma + 2.0	# something impossible
        two_sigma_m = (2 * sigma1 + sigma)/2
        # Iterate the following three equations 
        #  until there is no significant change in sigma 
        delsig = B*sin(sigma)*(cos(2*two_sigma_m)+0.25*B*(cos(sigma)*(2*cos(2*two_sigma_m)**2-1)-B/6*cos(2*two_sigma_m)*(-3+4*sin(sigma)**2)*(-3+4*cos(2*two_sigma_m)**2)))
        # two_sigma_m , delta_sigma
        count = 0
        while ( abs( two_sigma_m) > 0.01 and count<200) :
                two_sigma_m = 2 * sigma1 + sigma

                delta_sigma = B * sin(sigma) * ( cos(two_sigma_m) 
                        + (B/4) * (cos(sigma) * 
                        (-1 + 2 * pow( cos(two_sigma_m), 2 ) -  
                        (B/6) * cos(two_sigma_m) * \
                        (-3 + 4 * pow(sin(sigma), 2 )) *  \
                        (-3 + 4 * pow( cos (two_sigma_m), 2 ))))) 
                
                last_sigma = sigma
                sigma = (s / (b * A)) + delta_sigma
                count = count+1

        latitude2 = atan2 ( (sin(U1) * cos(sigma) + cos(U1) * sin(sigma) * cos(alpha1To2) ), \
                ((1-f) * sqrt( pow(Sinalpha, 2) +  \
                pow(sin(U1) * sin(sigma) - cos(U1) * cos(sigma) * cos(alpha1To2), 2))))

        lembda = atan2( (sin(sigma) * sin(alpha1To2 )), (cos(U1) * cos(sigma) -  \
                sin(U1) *  sin(sigma) * cos(alpha1To2)))

        C = (f/16) * cosalpha_sq * (4 + f * (4 - 3 * cosalpha_sq ))

        omega = lembda - (1-C) * f * Sinalpha *  \
                (sigma + C * sin(sigma) * (cos(two_sigma_m) + \
                C * cos(sigma) * (-1 + 2 * pow(cos(two_sigma_m),2) )))

        longitude2 = longitude1 + omega

        alpha21 = atan2 ( Sinalpha, (-sin(U1) * sin(sigma) +  \
                cos(U1) * cos(sigma) * cos(alpha1To2)))

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
        return longitude2, latitude2, alpha21 