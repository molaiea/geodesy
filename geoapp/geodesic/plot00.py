import datetime
import glob
import logging
import os
from numpy import *
from math import radians
import plotly.graph_objs as go
from plotly.offline import plot

def geodetic_to_geocentric(ellipsoid, latitude, longitude, height):

    φ = latitude
    λ = longitude
    sin_φ = sin(φ)
    a, rf = ellipsoid           
    e2 = 1 - (1 - 1 / rf) ** 2 
    n = a / sqrt(1 - e2 * sin_φ ** 2)
    r = (n + height) * cos(φ) 
    x = r * cos(λ)
    y = r * sin(λ)
    z = (n * (1 - e2) + height) * sin_φ
    return x, y, z

def plot3d(a, b, arr):
    phi = linspace(0, 2*pi,100)
    theta = linspace(0, pi,100)
    phi, theta=meshgrid(phi, theta)


    fig = go.Figure(go.Surface(
    x = cos(phi)*cos(theta)*7,
    y = cos(phi)*sin(theta)*7,
    z = 5 *sin(phi)))


    layout = go.Layout(
        title='Visualisation de la géodésique',
        autosize=False,
        width=800,
        height=800
        
        
    )
    fig = go.Figure(go.Surface(
    x = cos(phi)*cos(theta)*a,
    y = cos(phi)*sin(theta)*a,
    z = b *sin(phi)), layout=layout)

    spheroid = a,a/(a-b)
    for i in range(400):
        arr[i].x, arr[i].y, arr[i].z = geodetic_to_geocentric(spheroid, arr[i].Lat, arr[i].Long, 0)

    fig.add_scatter3d(
    x=[arr[i].x for i in range(400)],
    y=[arr[i].y for i in range(400)],
    z=[arr[i].z for i in range(400)],
    mode="lines"
)
    plot_div = plot(fig, output_type='div', include_plotlyjs=False)
    return plot_div