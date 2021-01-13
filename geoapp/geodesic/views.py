from django.shortcuts import render
from django.http import HttpResponse
from django.http import JsonResponse
import io
import urllib, base64
 
from math import *
import numpy as np
from . import plot00, directp
from . import problemedirect, problemeinverse
import os

from .forms import finalform, inverseform

# Create your views here.

def plott(request):
    print(request.POST)


def inverse(request):
    if request.method == 'POST':
        form = inverseform(request.POST)
        action = request.POST['action']
        if action == "Calculer":
            latitude, longitude, latitude0, longitude0 = 0, 0, 0, 0
            if form.is_valid():
                ellipsoid = form.cleaned_data.get("ellipsoid")
                a = form.cleaned_data.get("grand")
                b = form.cleaned_data.get("petit")
                if a==0 or not a or b==0 or not b:
                    if ellipsoid == "wgs":
                        a,b = 6378137,6356752.3142
                    elif ellipsoid == "grs":
                        a,b = 6378137,6356752.3141
                    elif ellipsoid == "clarke":
                        a,b =  	6378249.145,6356514.870
                latitude= form.cleaned_data.get("latitude")
                longitude= form.cleaned_data.get("longitude")
                latitude0= form.cleaned_data.get("latitude0")
                longitude0= form.cleaned_data.get("longitude0")
                
                s, az1, az2  = problemeinverse.inversefunction(a, b, latitude*pi/180, longitude*pi/180, latitude0*pi/180, longitude0*pi/180)
                azdirect = str(az1)+"°"
                azinverse = str(az2)+"°"
                distance = str(s)+"m"
                return render(request, 'inverse.html', {'form': form, 'az1': azdirect, 'az2': azinverse, 'distance': distance})
        elif action =="Visualiser":
            latitude, longitude, latitude0, longitude0 = 0, 0, 0, 0
            if form.is_valid():
                latitude= form.cleaned_data.get("latitude")
                latitude = latitude*pi/180
                longitude= form.cleaned_data.get("longitude")
                longitude = longitude*pi/180
                latitude0= form.cleaned_data.get("latitude0")
                latitude0 = latitude*pi/180
                longitude0= form.cleaned_data.get("longitude0")
                longitude0 = longitude*pi/180
                ellipsoid = form.cleaned_data.get("ellipsoid")
                a = form.cleaned_data.get("grand")
                b = form.cleaned_data.get("petit")
                if a==0 or b==0:
                    if ellipsoid == "wgs":
                        a,b = 6378137,6356752.3142
                    elif ellipsoid == "grs":
                        a,b = 6378137,6356752.3141
                    elif ellipsoid == "clarke":
                        a,b =  	6378249.145,6356514.870
                az1, az2, s  = problemeinverse.inversefunction(a, b, latitude, longitude, latitude0, longitude0)
                
                arr = problemedirect.geodesicpoints(a, b, latitude, longitude, az1*pi/180, s)
                return render(request, 'inverse.html', {'form': form, 'plot':plot00.plot3d(a, b, arr)})
    else:
        form = inverseform()
        return render(request, 'inverse.html', {'form': form})
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
    return atan(tanphi2)*180/pi, lamf*180/pi, alpha2*180/pi

def home(request):
    return render(request, 'home.html')
    


def direct(request):

    
    if request.method == 'POST':
        action = request.POST['action']
        if action =="Calculer":
        
            form = finalform(request.POST)
            latitude, longitude = 0, 0 
            if form.is_valid():
                latitude= form.cleaned_data.get("latitude")
                latitude = latitude*pi/180
                longitude= form.cleaned_data.get("longitude")
                longitude = longitude*pi/180
                ellipsoid = form.cleaned_data.get("ellipsoid")
                a = form.cleaned_data.get("grand")
                b = form.cleaned_data.get("petit")
                if a==0 or b==0:
                    if ellipsoid == "wgs":
                        a,b = 6378137,6356752.3142
                    elif ellipsoid == "grs":
                        a,b = 6378137,6356752.3141
                    elif ellipsoid == "clarke":
                        a,b =  	6378249.145,6356514.870
                azimut = form.cleaned_data.get("azimut")
                azimut = azimut*pi/180
                s = form.cleaned_data.get("distance_geodesique")

                lonf, latf, azf = directp.direct(a, b,latitude, longitude, azimut, s)
                if lonf<0:
                    lon2 = str(round(-lonf))+"° O"
                else:
                    lon2 = str(round(lonf))+"° E"
                if latf<0:
                    lat2 = str(round(latf))+"° S"
                else:
                    lat2 =str(round(latf))+"° N"
                azinverse = str(azf)+"°"
                #arr = [0 for i in range(400)]
                #print(Plot3DView.as_view())
                return render(request, 'direct.html', {'form': form, 'latitude': lat2, 'longitude': lon2, 'azimut': azinverse})
            
        elif action =="Visualiser":
            form = finalform(request.POST)
            latitude, longitude = 0, 0 
            if form.is_valid():
                latitude= form.cleaned_data.get("latitude")
                latitude = latitude*pi/180
                longitude= form.cleaned_data.get("longitude")
                longitude = longitude*pi/180
                ellipsoid = form.cleaned_data.get("ellipsoid")
                a = form.cleaned_data.get("grand")
                b = form.cleaned_data.get("petit")
                if a==0 or b==0:
                    if ellipsoid == "wgs":
                        a,b = 6378137,6356752.3142
                    elif ellipsoid == "grs":
                        a,b = 6378137,6356752.3141
                    elif ellipsoid == "clarke":
                        a,b =  	6378249.145,6356514.870
                azimut = form.cleaned_data.get("azimut")
                azimut = azimut*pi/180
                s = form.cleaned_data.get("distance_geodesique")
                
                #arr = [0 for i in range(400)]
                arr = problemedirect.geodesicpoints(a, b, latitude, longitude, azimut, s)
                print(arr)
                #print(Plot3DView.as_view())
                return render(request, 'direct.html', {'form': form, 'plot':plot00.plot3d(a, b, arr)})
    
    else:
        form = finalform()
        return render(request, 'direct.html', {'form': form})
