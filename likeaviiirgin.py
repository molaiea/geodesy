import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
import io
import matplotlib.pyplot as plt
from math import atan, sqrt, acos, asin

# introducin Point class (bash njmeo fiha les propriétés dstructure bash nshlu lkhdma)
class Point:
    Lat, Long, Az, Beta1, Beta0, W2, A1, B1, Sigma, alphaE, Sigma1, Sm, Dsigma, Sigmaf, x, y, z = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

#intro dlparameters dl ellipsoid
a = 6378137
b = 6356752.314245
f = (a-b)/a
#1ere eccentricitee e
e2 = sqrt((a**2-b**2)/a**2)
#2eme eccentricitee e'
e1 = sqrt((a**2-b**2)/b**2)

#coordonnees du point de depart lat long
phi1 = 2*pi/9
lambda1 = pi/6

#distance geodesique
s = 922300
#azimut du pt de depart
az1 = 0.2*pi
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

    Pi.Beta1 = atan((1-f)*tan(Pi.Lat))
    Pi.Beta0 = acos(cos(Pi.Beta1)*sin(az1))
    Pi.W2 = e1**2*(sin(Pi.Beta1)**2)
    Pi.Sigma1 = atan(tan(Pi.Beta1)/cos(az1))
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
    p.append(Pi)

#intro dyal sm spheroids ghi pr le test hna khdmt b WGS84
GRS80 = 6378137, 298.257222100882711
WGS84 = 6378137, 298.257223563

#function to convert geodetic coordinates(alpha, lambda, height) to geocentric(x, y, z)
def geodetic_to_geocentric(ellipsoid, latitude, longitude, height):
    """Return geocentric (Cartesian) Coordinates x, y, z corresponding to
    the geodetic coordinates given by latitude and longitude (in
    degrees) and height above ellipsoid. The ellipsoid must be
    specified by a pair (semi-major axis, reciprocal flattening).

    """
    φ = latitude
    λ = longitude
    sin_φ = sin(φ)
    a, rf = ellipsoid           # semi-major axis, reciprocal flattening
    e2 = 1 - (1 - 1 / rf) ** 2  # eccentricity squared
    n = a / sqrt(1 - e2 * sin_φ ** 2) # prime vertical radius
    r = (n + height) * cos(φ)   # perpendicular distance from z axis
    x = r * cos(λ)
    y = r * sin(λ)
    z = (n * (1 - e2) + height) * sin_φ
    return x, y, z

#loopin thro d points bash kula point nlqaw lih le x, y et z
for i in range(n):
    p[i].x, p[i].y, p[i].z = geodetic_to_geocentric(WGS84, p[i].Lat, p[i].Long, 0)

#ellipsoid plot
mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure(figsize = (10,10))
ax = fig.gca(projection='3d')
coefs = (a, a, b)  # Coefficients in a0/c x**2 + a1/c y**2 + a2/c z**2 = 1 
# Radii corresponding to the coefficients:
rx, ry, rz = coefs
u = linspace(0, 2 * pi, 100)
v = linspace(0, pi, 100)

# Cartesian coordinates that correspond to the spherical angles:
# (this is the equation of an ellipsoid):
x = rx * outer(cos(u), sin(v))
y = ry * outer(sin(u), sin(v))
z = rz * outer(ones_like(u), cos(v))

# Plot:
ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='g', alpha=0.5)
ax.plot([p[i].x for i in range(n)],
    [p[i].y for i in range(n)],
    [p[i].z for i in range(n)])

#coordonnees geocetriques dle pt de départ
xt, yt, zt =geodetic_to_geocentric(WGS84, phi1, lambda1, 0)
#phi lambda du pt darrivee
phif = 46.51404856*pi/180
lamf = 80.05730223*pi/180
#coordonnes geocentrique du pt darrivee
xf, yf, zf =geodetic_to_geocentric(WGS84, phif, lamf, 0)
#plotting pt de depart o darrivee
ax.scatter(xt, yt, zt, color='r')
ax.scatter(xf, yf, zf, color='r')
ax.legend()
plt.show()

#problem : geodesique mktmshish llpt darrivee