from math import sqrt, sin, cos, pi, atan, atan2, asin, acos, tan

def inversefunction(a, b, latitude1, longitude1, latitude2, longitude2):
      f = (a-b)/a		

      if (abs( latitude2 - latitude1 ) < 1e-8) and ( abs( longitude2 - longitude1) < 1e-8 ) :
            return 0.0, 0.0, 0.0


      b = a * (1.0 - f)

      beta1 = atan((1-f) * tan( latitude1 ))
      beta2 = atan((1-f) * tan( latitude2 ))

      delta_lambda = longitude2 - longitude1
      #variation_lambda = -4000000.0	#an impossibe value
      delta_u = delta_lambda

      # Iterate the following equations, 
      #  until there is no significant change in delta_lambda 
      sqr_sin_sigma = pow( cos(beta2) * sin(delta_lambda), 2) + pow( (cos(beta1) * sin(beta2) - sin(beta1) *  cos(beta2) * cos(delta_lambda) ), 2 )
      Sin_sigma = sqrt( sqr_sin_sigma )
      Cos_sigma = sin(beta1) * sin(beta2) + cos(beta1) * cos(beta2) * cos(delta_lambda)
      sigma = atan2( Sin_sigma, Cos_sigma )
      Sin_alpha = cos(beta1) * cos(beta2) * sin(delta_lambda) / sin(sigma)
      alpha = asin( Sin_alpha )
      Cos2sigma_m = cos(sigma) - (2 * sin(beta1) * sin(beta2) / pow(cos(alpha), 2) )
      C = (f/16) * pow(cos(alpha), 2) * (4 + f * (4 - 3 * pow(cos(alpha), 2)))

      variation_lambda = delta_lambda

      delta_lambda = delta_u + (1-C) * f * sin(alpha) * (sigma + C * sin(sigma) * \
                  (Cos2sigma_m + C * cos(sigma) * (-1 + 2 * pow(Cos2sigma_m, 2) )))

      while variation_lambda - delta_lambda > 1.0e-9  :

            sqr_sin_sigma = pow( cos(beta2) * sin(delta_lambda), 2) + \
                  pow( (cos(beta1) * sin(beta2) - \
                  sin(beta1) *  cos(beta2) * cos(delta_lambda) ), 2 )

            Sin_sigma = sqrt( sqr_sin_sigma )

            Cos_sigma = sin(beta1) * sin(beta2) + cos(beta1) * cos(beta2) * cos(delta_lambda)

            sigma = atan2( Sin_sigma, Cos_sigma )

            Sin_alpha = cos(beta1) * cos(beta2) * sin(delta_lambda) / sin(sigma)
            alpha = asin( Sin_alpha )

            Cos2sigma_m = cos(sigma) - (2 * sin(beta1) * sin(beta2) / pow(cos(alpha), 2) )

            C = (f/16) * pow(cos(alpha), 2) * (4 + f * (4 - 3 * pow(cos(alpha), 2)))

            variation_lambda = delta_lambda

            delta_lambda = delta_u + (1-C) * f * sin(alpha) * (sigma + C * sin(sigma) * \
                  (Cos2sigma_m + C * cos(sigma) * (-1 + 2 * pow(Cos2sigma_m, 2) )))

      w_carre = pow(cos(alpha),2) * (a*a-b*b) / (b*b)

      A = 1 + (w_carre/16384) * (4096 + w_carre * (-768 + w_carre * (320 - 175 * w_carre)))

      B = (w_carre/1024) * (256 + w_carre * (-128+ w_carre * (74 - 47 * w_carre)))

      delta_sigma = B * Sin_sigma * (Cos2sigma_m + (B/4) * \
            (Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2) ) - \
            (B/6) * Cos2sigma_m * (-3 + 4 * sqr_sin_sigma) * \
            (-3 + 4 * pow(Cos2sigma_m,2 ) )))

      s = b * A * (sigma - delta_sigma)

      alpha12 = atan2( (cos(beta2) * sin(delta_lambda)), \
            (cos(beta1) * sin(beta2) - sin(beta1) * cos(beta2) * cos(delta_lambda)))

      alpha21 = atan2( (cos(beta1) * sin(delta_lambda)), \
            (-sin(beta1) * cos(beta2) + cos(beta1) * sin(beta2) * cos(delta_lambda)))

      if ( alpha12 < 0.0 ) : 
            alpha12 =  alpha12 + 2*pi
      if ( alpha12 > 2*pi ) : 
            alpha12 = alpha12 - 2*pi

      alpha21 = alpha21 + 2*pi / 2.0
      if ( alpha21 < 0.0 ) : 
            alpha21 = alpha21 + 2*pi
      if ( alpha21 > 2*pi ) : 
            alpha21 = alpha21 - 2*pi

      alpha12 = alpha12*180/pi
      alpha21 = alpha21*180/pi
      return round(s), round(alpha12),  round(alpha21)
