from django.shortcuts import render
from django.http import HttpResponse
from django.http import JsonResponse
import io
import urllib, base64
 
from math import *
import numpy as np
from . import plot00
from . import problemedirect
import os
import plotly.graph_objs as go
from plotly.offline import plot
from django.views.generic import TemplateView
from django.shortcuts import render, redirect
from .forms import finalform, LatitudeField

# Create your views here.

def plott(request):
    print(request.POST)



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

    # alpha2 = atan(sin(alphae)/(cos(beta1)*cos(sigma)*cos(alpha1)-sin(beta1)*sin(sigma)))
    return atan(tanphi2), lamf

def home(request):
    return render(request, 'home.html')
    
def inverse(request):
    return render(request, 'inverse.html')

def direct(request):

    
    print(request.method)
    if request.method == 'POST':
        submit= request.POST.get("submit")
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
            latf, lonf = geodesic(longitude, latitude, azimut, s, a, b)
            #arr = [0 for i in range(400)]
            arr = problemedirect.geodesicpoints(a, b, latitude, longitude, azimut, s)
            
            #print(Plot3DView.as_view())
        return render(request, 'direct.html', {'form': form, 'latitude': latf*180/pi, 'longitude': lonf*180/pi, 'plot':plot00.plot3d(a, b, arr)})
    else:
        form = finalform()
    return render(request, 'direct.html', {'form': form})
