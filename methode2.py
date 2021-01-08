# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 10:33:55 2020

@author: Moussaoui
"""
#test youssef
import math
#Paramètres de l'ellipsoïde WGS 84
ac= 6378137.0 
bc= 6356752.314245 
f= (ac - bc) / ac 
e2 = math.sqrt((ac**2 - bc**2) / bc**2) 
Phi1=int(input('entrer latitde du point de depart'))
Phi2=int(input('entrer latitude du point darrivée'))
Lambda1=int(input('entrer longitude du point de départ'))
Lambda2=int(input('entrer longitude du point darrivée'))
#Problème inverse
#onversion de la longitude et la latitude du point de départ en radian
Pi = math.pi
 #Affichage de la conversion en degré
#♥Conversion des angles en radian pour les calculs
Phi1 = Phi1 * (Pi / 180) 
Phi2 = Phi2 * (Pi / 180)
Lambda1 = Lambda1*(Pi / 180)
Lambda2 = Lambda2*(Pi / 180)
#Calcul des latitudes réduites au point du départ
B1 = math.atan((1 - f) * math.tan(Phi1))
B2 = math.atan((1 - f) * math.tan(Phi2))
#Calcul de la différence de longitude sur l'ellipsoïde
Delta1 = Lambda2 - Lambda1
#Calcul de la différence de longitude sur la sphère
#Calcul des grandeurs suivantes
Sigma = math.asin(math.sqrt((math.cos(B2) * math.sin(Delta1))**2 + (math.cos(B1) * math.sin(B2) - math.sin(B1) * math.cos(B2) * math.cos(Delta1))**2))
AE = math.asin((math.cos(B1) * math.cos(B2) * math.sin(Delta1)) / math.sin(Sigma))
Sm = math.acos(math.cos(Sigma) - (2 * math.sin(B1) * math.sin(B2) / (math.cos(AE)**2))) / 2
#Calcul de la constante de Vicenty
C = (f / 16) * (math.cos(AE)**2) * (4 + f * (4 - 3 * (math.cos(AE)**2)))
#Calcul de la diférence de longitude sur la sphère auxiliaire
Delta2 = Delta1 + (1 - C) * f * math.sin(AE) * (Sigma + C * math.sin(Sigma) * (math.cos(2 * Sm) + C * math.cos(Sigma) * (-1 + 2 * (math.cos(2 * Sm)**2))))
#Itérations jusqu'a que la valeur de la différence de longitude sur la sphère soit négligeable
Deltau= Delta1
while abs(Deltau) < 0.000001: 
                Sigma = math.asin(math.sqrt((math.Cos(B2) * math.sin(Delta2))**2 + (math.cos(B1) * math.sin(B2) - math.sin(B1) * math.cos(B2) * math.cos(Delta2))**2))
                AE = math.asin((math.cos(B1) * math.cos(B2) * math.sin(Delta2)) / math.sin(Sigma))
                Sm = math.acos(math.cos(Sigma) - (2 * math.sin(B1) * math.sin(B2) / (math.cos(AE)**2))) / 2
                C = (f / 16) * (math.cos(AE)**2) * (4 + f * (4 - 3 * (math.cos(AE)**2)))
                Deltau = Delta1 + (1 - C) * f * math.sin(AE) * (Sigma + C * math.sin(Sigma) * (math.cos(2 * Sm) + C * math.cos(Sigma) * (-1 + 2 * (math.cos(2 * Sm)**2))))
                Delta1=Deltau
print (Deltau)
