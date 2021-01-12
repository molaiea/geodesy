import math
import numpy as np

def calculateGeographicalPositionFromRangeBearing(a, b, latitude1, longitude1, alpha1To2, s ) :


        f = (a-b)/a			# metres

        piD4 = math.atan( 1.0 )
        two_pi = piD4 * 8.0

        latitude1    = latitude1    * piD4 / 45.0
        longitude1 = longitude1 * piD4 / 45.0
        alpha1To2 = alpha1To2 * piD4 / 45.0
        if ( alpha1To2 < 0.0 ) : 
                alpha1To2 = alpha1To2 + two_pi
        if ( alpha1To2 > two_pi ) : 
                alpha1To2 = alpha1To2 - two_pi

        TanU1 = (1-f) * math.tan(latitude1)
        U1 = math.atan( TanU1 )
        sigma1 = math.atan2( TanU1, math.cos(alpha1To2) )
        Sinalpha = math.cos(U1) * math.sin(alpha1To2)
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
        delsig = B*math.sin(sigma)*(math.cos(2*two_sigma_m)+0.25*B*(math.cos(sigma)*(2*math.cos(2*two_sigma_m)**2-1)-B/6*math.cos(2*two_sigma_m)*(-3+4*math.sin(sigma)**2)*(-3+4*math.cos(2*two_sigma_m)**2)))
        # two_sigma_m , delta_sigma
        count = 0
        while ( abs( two_sigma_m) > 0.01 and count<200) :
                two_sigma_m = 2 * sigma1 + sigma

                delta_sigma = B * math.sin(sigma) * ( math.cos(two_sigma_m) 
                        + (B/4) * (math.cos(sigma) * 
                        (-1 + 2 * math.pow( math.cos(two_sigma_m), 2 ) -  
                        (B/6) * math.cos(two_sigma_m) * \
                        (-3 + 4 * math.pow(math.sin(sigma), 2 )) *  \
                        (-3 + 4 * math.pow( math.cos (two_sigma_m), 2 ))))) 
                
                last_sigma = sigma
                sigma = (s / (b * A)) + delta_sigma
                count = count+1

        latitude2 = math.atan2 ( (math.sin(U1) * math.cos(sigma) + math.cos(U1) * math.sin(sigma) * math.cos(alpha1To2) ), \
                ((1-f) * math.sqrt( math.pow(Sinalpha, 2) +  \
                pow(math.sin(U1) * math.sin(sigma) - math.cos(U1) * math.cos(sigma) * math.cos(alpha1To2), 2))))

        lembda = math.atan2( (math.sin(sigma) * math.sin(alpha1To2 )), (math.cos(U1) * math.cos(sigma) -  \
                math.sin(U1) *  math.sin(sigma) * math.cos(alpha1To2)))

        C = (f/16) * cosalpha_sq * (4 + f * (4 - 3 * cosalpha_sq ))

        omega = lembda - (1-C) * f * Sinalpha *  \
                (sigma + C * math.sin(sigma) * (math.cos(two_sigma_m) + \
                C * math.cos(sigma) * (-1 + 2 * math.pow(math.cos(two_sigma_m),2) )))

        longitude2 = longitude1 + omega

        alpha21 = math.atan2 ( Sinalpha, (-math.sin(U1) * math.sin(sigma) +  \
                math.cos(U1) * math.cos(sigma) * math.cos(alpha1To2)))

        alpha21 = alpha21 + two_pi / 2.0
        if ( alpha21 < 0.0 ) :
                alpha21 = alpha21 + two_pi
        if ( alpha21 > two_pi ) :
                alpha21 = alpha21 - two_pi

        latitude2       = latitude2       * 45.0 / piD4
        longitude2    = longitude2    * 45.0 / piD4
        alpha21    = alpha21    * 45.0 / piD4

        return longitude2, latitude2, alpha21 