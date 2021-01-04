import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats


# Table of contents #

# (0) ##### Simplifications ##### #
# (0.1) vector3(x,y,z) ... vector 3D simplification
# (0.2) vector4(x,y,z) ... vector 4D simplification (e.g. Quaternion)

# (1) ##### simple functions ##### #
# (1.1) radiant to degree conversion
# (1.2) degree to radiant conversion


# (2) ##### statistic functions ##### #
# (2.1) plot of normal distribution

# (3) ##### rotations, quaternion, euler #####
# (3.1) Rotmat2Euler( rotmat ) --- Rotation Matrix to Euler angles
# (3.2) vec2skewmat( x ): --- vector to skew symmetric martix
# (3.3) quat_mult( p, q ): --- quaterion multiplication --> concatenation of rotations
# (3.4) quat_change( sigK ): --- change in rotation as quaternion
# (3.5) rpy = Quat2Euler( q ): --- Quaternion to Euler angles 
# (3.6) q = Euler2Quat(r,p,y): --- Euler angles to Quaternion
# (3.7) R = quat2Rotmat( q ) --- Quaternion to Rotation matrix conversion
# (3.8) q = Rotmat2quat( R ) --- Rotation Matrix to Quaternion conversion

# (4) ##### transformations ##### #
# (4.1) ECEF-Frame to LLA Transformation, different Ellipsoids avaivable
# (4.2) Reference Ellipsoid function
# (4.3) ell2utm(lat,lon) ... Ellipsoidical Coordinates to UTM Transformation
# (4.4) utm2ell(N,E,Zone) ... UTM to LLA conversion
# (4.5) ell2xyz(lat,lon,h): --- lat, lon, h to ECEF (XYZ) coordinates

# (5) Gravitiy Models
# (5.1) Gravity Model ECEF - Frame


# (0) ##### Simplifications ##### #
def vector2(x,y):
    return np.array((x,y), dtype=float)

def vector3(x,y,z):
    return np.array((x,y,z), dtype=float)

# (0.2) 4D float vector
def vector4(w,x,y,z):
    return np.array((w,x,y,z), dtype=float)

# (0.2) 4D float vector
def vector5(x1,x2,x3,x4,x5):
    return np.array((x1,x2,x3,x4,x5), dtype=float)

# (0.3) 4D float vector
def vector6(x1,x2,x3,x4,x5,x6):
    return np.array((x1,x2,x3,x4,x5,x6), dtype=float)  


# (1) ##### Simple Functions ##### #

# (1.1) radiant to degree conversion
def rad2deg(angle_rad):
# Input: angle in radiant
# Output: angle in degree
    return ( angle_rad * (180/np.pi) )

# (1.2) degree to radiant conversion
def deg2rad(angle_deg):
# Input: angle in degree
# Output: angle in radiant
    return ( angle_deg * (np.pi/180) )


# (2) ##### statistic functions ##### #

# (2.1) plot of normal distribution
def plot_normal_distr_1D(mu, sigma):

    x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
    y = stats.norm.pdf(x, mu, sigma) 

    # scale to maximum at 1
    y = y/y.max()

    plt.plot(x, y)

# (3) ##### rotations, quaternion, euler #####

# (3.1) Rotmat2Euler( rotmat ) --- Rotation Matrix to Euler angles
def Rotmat2Euler( rotmat ):
# ############################################################################
# Function computes the Euler Angles from given rotation matrix
# ----------------------------------------------------------------------------
# Input:
# rotmat (double 3x3)...matrix
#-----------------------------------------------------------------------------
# @author: Felix Esser
# @date: 14.07.2020
# @mail: s7feesse@uni-bonn.de
# @ literature: Förstner and Wrobel (2016), Photogrammetric Computer Vision
# ############################################################################


    roll  = np.arctan2( rotmat[2,1], rotmat[2,2] )
    pitch = np.arctan2( -rotmat[2,0], np.sqrt(rotmat[2,1]**2 + rotmat[2,2]**2) )
    yaw   = np.arctan2( rotmat[1,0], rotmat[0,0] )

    rpy = vector3(roll, pitch, yaw)

    return rpy

# (3.2) vec2skewmat( x ): --- vector to skew symmetric martix
def vec2skewmat( x ):
# ############################################################################
# Function computes the skew matrix from a given vector
# ----------------------------------------------------------------------------
# Input:
# v (double 3x1)...vector
# Output:
# S_x (double 3x3)...matrix
#-----------------------------------------------------------------------------
# @author: Felix Esser
# @date: 15.07.2020
# @mail: s7feesse@uni-bonn.de
# @literature: Förstner & Wrobel (2016), p.201
# ############################################################################

    S_x = np.zeros((3,3))

    S_x[0,0] =    0
    S_x[0,1] = -x[2]
    S_x[0,2] =  x[1]

    S_x[1,0] =  x[2]
    S_x[1,1] =    0
    S_x[1,2] = -x[0]

    S_x[2,0] = -x[1]
    S_x[2,1] =  x[0]
    S_x[2,2] =    0

    return S_x

# (3.3) quat_mult( p, q ): --- quaterion multiplication 
def quat_mult( q, r ):
# ############################################################################
# function computes the product of two quaternions, which applys the change in Orientation
# see Groves 2013 Appendix_E.pdf p.10
# author: Felix Esser
# date: 14.07.20
# multiplication: p = q * r, update of orientation with quaternions
#
# literatur: Förstner & Wrobel 2016, Photogrammetric Computer Vision, p.334, formula  (8.43)
# ############################################################################
    
    p = np.zeros((4))
    p[0] = q[0]*r[0] - q[1]*r[1] - q[2]*r[2] - q[3]*r[3]
    p[1] = q[1]*r[0] + q[0]*r[1] - q[3]*r[2] + q[2]*r[3]
    p[2] = q[2]*r[0] + q[3]*r[1] + q[0]*r[2] - q[1]*r[3]
    p[3] = q[3]*r[0] - q[2]*r[1] + q[1]*r[2] + q[0]*r[3]

    return p


# (3.4) quat_change( sigK ): --- change in rotation as quaternion
def quat_change( sigK ):
# ############################################################################
#  computes the Quaternion that represents the change in rotation
#  
#  date: 14.07.20
#  author: Felix Esser
#  E-mail: s7feesse@uni-bonn.de
#  literature: Grooves 2013, Appendix
# ############################################################################
    n = np.linalg.norm( sigK )

    rK = np.zeros((4))

    if (n > 10**-10):
        rK[0] =  np.cos(n/2)
        rK[1:4]  = (sigK/n) * (np.sin(n/2))
    else:
        rK[0] = 1 - (1/(2*(2**2)))*(n**2) + (1/(24*(2**4)))*(n**4) - (1/(720*(2**6)))*(n**6)
        rK[1:4] = sigK * (0.5 - (1/(6*(2**3)))*(n**2) + (1/(120*(2**5)))*(n**4) - (1/(5040*(2**7)))*(n**6) )
    return rK

# matlab

# n=norm(dsig); % norm
# if n>1e-10
#    r=[cos(n/2) ; dsig ./n.* sin(n/2)]; %(3.116)
# else
#    den=@(x) factorial(x)*2^x;   
#    r=[1-1/den(2)*n.^2 + 1/den(4)*n.^4 - 1/den(6)*n.^6;
#        dsig.*repmat(0.5 - 1/den(3)*n.^2 + 1/den(5)*n.^4 - 1/den(7)*n.^6,3,1)]'; %(3.117)
# end

# (3.5) Quat2Euler( q ): --- Quaternion to Euler angles 
def Quat2Euler( q ):
# ############################################################################
#  computes roll pitch yaw angle from quaternion 
# 
#  date: 13.07.20
#  author: Felix Esser
#  E-mail: s7feesse@uni-bonn.de
# ############################################################################

    q0 = q[0]
    q1 = q[1]
    q2 = q[2]
    q3 = q[3]

    roll = np.arctan2(2*(q0*q1+q2*q3),(1-2*q1**2-2*q2**2))
    pitch = np.arcsin(2*(q0*q2-q1*q3))
    yaw = np.arctan2(2*(q0*q3+q1*q2),(1-2*q2**2-2*q3**2))

    return roll, pitch, yaw

# (3.6) Euler2Quat(r,p,y): --- Euler angles to Quaternion
def Euler2Quat(r,p,y):
# ############################################################################
# 
# Function computes the Quaternion from Euler angles 
#  
# INPUT:    r,p,y:    orientation: roll, pitch und yaw ( EULER Angles )
# 
# OUTPUT:   q [4x1]:  orientation: quaternion
# 
# author: Felix Esser
# date:   13.07.20
# 
# ############################################################################

    # % Aufstellen des Quaternion
    q = np.zeros( ( 4, 1 ) )

    q[0,0] = np.cos(r/2)*np.cos(p/2)*np.cos(y/2) + np.sin(r/2)*np.sin(p/2)*np.sin(y/2)
    q[1,0] = np.sin(r/2)*np.cos(p/2)*np.cos(y/2) - np.cos(r/2)*np.sin(p/2)*np.sin(y/2)
    q[2,0] = np.cos(r/2)*np.sin(p/2)*np.cos(y/2) + np.sin(r/2)*np.cos(p/2)*np.sin(y/2)
    q[3,0] = np.cos(r/2)*np.cos(p/2)*np.sin(y/2) - np.sin(r/2)*np.sin(p/2)*np.cos(y/2)

    # old matlab code:
    # q = [ np.cos(r/2)*np.cos(p/2)*np.cos(y/2) + np.sin(r/2)*np.sin(p/2)*np.sin(y/2)
    #       np.sin(r/2)*np.cos(p/2)*np.cos(y/2) - np.cos(r/2)*np.sin(p/2)*np.sin(y/2)
    #       np.cos(r/2)*np.sin(p/2)*np.cos(y/2) + np.sin(r/2)*np.cos(p/2)*np.sin(y/2)
    #       np.cos(r/2)*np.cos(p/2)*np.sin(y/2) - np.sin(r/2)*np.sin(p/2)*np.cos(y/2) ];

    return q

# (3.7) R = quat2Rotmat( q ) Quaternion to Rotation matrix conversion
def quat2Rotmat( q ):
# ############################################################################
# Determine the rotation matrix R that corresponds to Quaternion
# ----------------------------------------------------------------------------
# Input:
# q 
# Output:
# R
#-----------------------------------------------------------------------------
# @author: Felix Esser
# @date: 14.07.2020
# @mail: s7feesse@uni-bonn.de
# literature: Förstner & Wrobel (2016), Photogrammetric Computer Vision, p. 335, formula 8.55
# ############################################################################

    a = q[0]
    b = q[1]
    c = q[2]
    d = q[3]

    R = np.zeros((3,3))
          
    # DCM Matrix (body => navigation!!):    
    # DCM = n*[(a^2+b^2-c^2-d^2)   2*(b*c-a*d)         2*(b*d+a*c);...
    #          2*(b*c+a*d)        (a^2-b^2+c^2-d^2)   2*(c*d-a*b);...
    #           2*(b*d-a*c)         2*(c*d+a*b)       (a^2-b^2-c^2+d^2)];       


    R[0,0] = a**2+b**2-c**2-d**2
    R[0,1] = 2*(b*c-a*d)
    R[0,2] = 2*(b*d+a*c)

    R[1,0] = 2*(b*c+a*d)
    R[1,1] = (a**2-b**2+c**2-d**2)
    R[1,2] = 2*(c*d-a*b)

    R[2,0] = 2*(b*d-a*c) 
    R[2,1] = 2*(c*d+a*b)
    R[2,2] = (a**2-b**2-c**2+d**2)

    R = 1/(np.linalg.norm(q)**2)*R

    return R

    # matlab code 
    # DCM = [(a^2+b^2-c^2-d^2)   2*(b*c-a*d)         2*(b*d+a*c);...
    #           2*(b*c+a*d)        (a^2-b^2+c^2-d^2)   2*(c*d-a*b);...
    #           2*(b*d-a*c)         2*(c*d+a*b)       (a^2-b^2-c^2+d^2)];

# (3.8) q = Rotmat2quat( R ) Rotation Matrix to Quaternion conversion
def Rotmat2quat( R ):

    # to euler
    rpy = Rotmat2Euler( R )

    # euler to quat 
    q = Euler2Quat( rpy[0], rpy[1], rpy[2] )

    return q




# (4) ##### transformations ##### #

# (4.1) ECEF-Frame to LLA Transformation
def ecef2lla(x, y, z, *left):

    # ecef2lla - convert earth-centered earth-fixed (ECEF)
    #            cartesian coordinates to latitude, longitude,
    #            and altitude (above WGS84 ellipsoid)

    # INPUT 
    # x = ECEF X-coordinate (m)
    # y = ECEF Y-coordinate (m)
    # z = ECEF Z-coordinate (m)

    # OUTPUT
    # lat = geodetic latitude (radians) B
    # lon = longitude (radians) L
    # alt = height above WGS84 ellipsoid (m)

    # modified: 08.07.20 by Felix Esser
    # - vector transformation now possible

    a = 6378137
    e = 8.1819190842622e-2

    #% calculations
    b = np.sqrt(a**2*(1-e**2))
    ep = np.sqrt((a**2-b**2)/b**2)

    p = np.sqrt(np.power(x, 2) + np.power(y, 2))
    th = np.arctan2(a*z, b*p)
    lon = np.arctan2(y, x)

    # lattitude --> Breitengrad B
    lat = np.arctan2((z+np.power(ep,2)*b*np.power(np.sin(th),3)), (p-np.power(e,2)*a*np.power(np.cos(th),3)))


    N = a / np.sqrt( 1 - np.power(e,2) * np.power(np.sin(lat),2) )

    alt = p/np.cos(lat)-N

    # longitude --> Längengrad L
    lon = np.fmod(lon, 2*np.pi)

    #% correct for numerical instability in altitude near exact poles:
    #% (after this correction, error is about 2 millimeters, which is about
    #% the same as the numerical precision of the overall function)

    if ( type(x) == type(np.array([1])) ): # array of coordinates
        for i in range(0,len(x)):
            if (np.fabs(x[i]) < 1 and np.fabs(y[i]) < 1):
                alt[i] = np.fabs(z[i])-b
    else: # single coordinate
        if (np.fabs(x) < 1 and np.fabs(y) < 1):
                alt = np.fabs(z)-b

    return (lat, lon, alt)


# (4.2) Reference Ellipsoid function
def refell(sys):
# modified by: Felix Esser, 20.07.20
# REFELL  Computes reference ellispoid parameters.
#   TOPEX Reference: <http://topex-www.jpl.nasa.gov/aviso/text/general/news/hdbk311.htm#CH3.3>
# Version: 1 Jun 04
# Useage:  [a,b,e2,finv]=refell(type)
# Input:   type - reference ellipsoid type (char)
#                 CLK66 = Clarke 1866
#                 GRS67 = Geodetic Reference System 1967
#                 GRS80 = Geodetic Reference System 1980
#                 WGS72 = World Geodetic System 1972
#                 WGS84 = World Geodetic System 1984
#                 ATS77 = Quasi-earth centred ellipsoid for ATS77
#                 NAD27 = North American Datum 1927 (=CLK66)
#                 NAD83 = North American Datum 1927 (=GRS80)
#                 INTER = International
#                 KRASS = Krassovsky (USSR)
#                 MAIRY = Modified Airy (Ireland 1965/1975)
#                 TOPEX = TOPEX/POSEIDON ellipsoid
# Output:  a    - major semi-axis (m)
#          b    - minor semi-axis (m)
#          e2   - eccentricity squared
#          finv - inverse of flattening

    if (sys =='CLK66' or sys =='NAD27'):
      a=6378206.4;
      finv=294.9786982;
    elif sys =='GRS67':
      a=6378160.0;
      finv=298.247167427;
    elif (sys =='GRS80' or sys =='NAD83'):
      a=6378137.0;
      finv=298.257222101;
    elif (sys =='WGS72'):
      a=6378135.0;
      finv=298.26;
    elif (sys =='WGS84'):
      a=6378137.0;
      finv=298.257223563;
    elif sys =='ATS77':
      a=6378135.0;
      finv=298.257;
    elif sys =='KRASS':
      a=6378245.0;
      finv=298.3;
    elif sys =='INTER':
      a=6378388.0;
      finv=297.0;
    elif sys =='MAIRY':
      a=6377340.189;
      finv=299.3249646;
    elif sys =='TOPEX':
      a=6378136.3;
      finv=298.257;

    f=1/finv;
    b=a*(1-f);
    e2=1-(1-f)**2;

    return a, b, e2, finv


# (4.3) Ellipsoidical Coordinates to UTM Transformation
def ell2utm(lat,lon,Zone):
# function [N,E,Zone]=ell2utm(lat,lon,a,e2,lcm)
# ELL2UTM  Converts ellipsoidal coordinates to UTM.
#   UTM northing and easting coordinates in a 6 degree
#   system.  Zones begin with zone 1 at longitude 180E
#   to 186E and increase eastward.  Formulae from E.J.
#   Krakiwsky, "Conformal Map Projections in Geodesy",

# Useage:  [N,E,Zone]=ell2utm(lat,lon,a,e2,lcm)
#          [N,E,Zone]=ell2utm(lat,lon,a,e2)
#          [N,E,Zone]=ell2utm(lat,lon,lcm)
#          [N,E,Zone]=ell2utm(lat,lon)
# Input:   lat - vector of latitudes (rad)
#          lon - vector of longitudes (rad)
#          a   - ref. ellipsoid major semi-axis (m); default GRS80
#          e2  - ref. ellipsoid eccentricity squared; default GRS80
#          lcm - central meridian; default = standard UTM def'n
# Output:  N   - vector of UTM northings (m)
#          E   - vector of UTM eastings (m)
#          Zone- vector of UTM zones

    # variable to save UTM north east
    [a, b, e2, finv] = refell('GRS80')

    lcm = (Zone*6-183) * (np.pi/180)    # ML --> lcm=deg2rad(Zone*6-183);

    ko=0.9996                 # Scale factor
    No=np.zeros( lat.shape )  # False northing (north)
    # No(lat<0)=1e7           # False northing (south)
    np.where(lat < 0, No, 1e7)
    Eo=500000                 # False easting

    lam=lon-lcm

    np.where(lam >= np.pi, lam, lam - 2*np.pi) # lam=lam-(lam >= pi)*(2*pi)
      
    #fprintf('\nZones\n');
    #fprintf('%3d\n',Zone');
    #fprintf('\nCentral Meridians\n');
    #fprintf('%3d %2d %9.6f\n',rad2dms(lcm)');
    #fprintf('\nLongitudes wrt Central Meridian\n');
    #fprintf('%3d %2d %9.6f\n',rad2dms(lam)');

    f = 1 - np.sqrt(1 - e2)
    RN = a / ((1 - e2 * np.sin(lat)**2))**0.5
    RM = a*(1-e2)/(1-e2 * np.sin(lat)**2)**1.5
    t = np.tan(lat)
    h = np.sqrt((e2*np.cos(lat)**2) / (1 - e2))
    n = f/(2-f)

    a0=1+n**2 / 4+n**4 / 64
    a2=1.5*(n-n**3/8)
    a4=15/16*(n**2-n**4/4)
    a6=35/48*n**3
    a8=315/512*n**4

    s=a/(1+n)*(a0*lat-a2*np.sin(2*lat)+a4*np.sin(4*lat)-a6*np.sin(6*lat)+a8*np.sin(8*lat))

    E1=lam * np.cos(lat)
    E2=lam**3 * np.cos(lat)**3/6 * (1-t**2+h**2)
    E3=lam**5 * np.cos(lat)**5/120 * (5-18*t**2+t**4+14*h**2-58*t**2*h**2+13*h**4+4*h**6-64*t**2*h**4-24*t**2*h**6)
    E4=lam**7 * np.cos(lat)**7/5040*(61-479*t**2+179*t**4-t**6)

    E=Eo + ko*RN*(E1 + E2 + E3 + E4)

    N1=lam**2/2 * np.sin(lat) * np.cos(lat)
    N2=lam**4/24 * np.sin(lat) * np.cos(lat)**3*(5-t**2+9*h**2+4*h**4)
    N3=lam**6/720 * np.sin(lat) * np.cos(lat)**5*(61-58*t**2+t**4+270*h**2-330*t**2*h**2+445*h**4+324*h**6-680*t**2*h**4+88*h**8-600*t**2*h**6-192*t**2*h**8)
    N4=lam**8/40320*np.sin(lat)*np.cos(lat)**7*(1385-311*t**2+543*t**4-t**6)

    N=No + ko*RN*(s/RN + N1 + N2 + N3 + N4)

    return N,E


def utm2ell(N,E,Zone):
    # get parameters of grs80 Ellipsoid
    [a, b, e2, finv] = refell('GRS80')
    f = 1 / finv

    #----- Remove false northings and eastings
    No=np.zeros(N.shape)   # False northing (north)
    #No(Zone<0)=1e7      # False northing (south)
    Eo=500000           # False easting
    N=N-No
    E=E-Eo
    Zone=np.absolute(Zone)   # Remove negative zone indicator for southern hemisphere

    #----- Foot point latitude
    ko=0.9996           # UTM scale factor
    lat1=N/ko/a
    dlat=1

    lcm=deg2rad(np.absolute(Zone)*6-183)

    while np.amax(np.absolute(dlat))>1e-12:
        A0=1-(e2/4)-(e2**2*3/64)-(e2**3*5/256)-(e2**4*175/16384)
        A2=(3/8)*( e2+(e2**2/4)+(e2**3*15/128)-(e2**4*455/4096) )
        A4=(15/256)*( e2**2+(e2**3*3/4)-(e2**4*77/128) )
        A6=(35/3072)*( e2**3-(e2**4*41/32) )
        A8=-(315/131072)*e2**4
        f1=a*( A0*lat1-A2*np.sin(2*lat1)+A4*np.sin(4*lat1)-A6*np.sin(6*lat1)+A8*np.sin(8*lat1) )-N/ko
        f2=a*( A0-2*A2*np.cos(2*lat1)+4*A4*np.cos(4*lat1)-6*A6*np.cos(6*lat1)+8*A8*np.cos(8*lat1) )

        dlat=-f1/f2
        lat1=lat1+dlat

    RN=a/(1-e2*np.sin(lat1)**2)**0.5

    RM=a*(1-e2)/(1-e2*np.sin(lat1)**2)**1.5
    h2=e2*np.cos(lat1)**2/(1-e2)
    t=np.tan(lat1)


    E0=E/ko/RN
    E1=E0
    E2=E0**3/6*(1+2*t**2+h2)
    E3=E0**5/120*(5+6*h2+28*t**2-3*h2**2+8*t**2*h2+24*t**4-4*h2**3+4*t**2.*h2**2+24*t**2*h2**3)
    E4=E0**7/5040*(61 + 662*t**2 + 1320*t**4 + 720*t**6)
    lon=(1/np.cos(lat1))*(E1-E2+E3-E4)+lcm

    E0=E/ko
    N1=(t*E0**2)/(2*RM*RN)
    N2=(t*E0**4)/(24*RM*RN**3)*(5+3*t**2+h2-4*h2**2-9*h2*t**2)
    N3=(t*E0**6)/(720*RM*RN**5)*(61-90*t**2+46*h2+45*t**4-252*t**2*h2-5*h2**2+100*h2**3-66*t**2*h2**2-90*t**4.*h2+88*h2**4+225*t**4.*h2**2+84*t**2*h2**3 - 192*t**2*h2**4)
    N4=(t*E0**8)/(40320*RM*RN**7)*(1385+3633*t**2+4095*t**4+1575*t**6)
    lat=lat1-N1+N2-N3+N4

    return lat, lon


def ell2xyz(lat,lon,h):
    [a, b, e2, finv] = refell('GRS80')

    v=a/np.sqrt(1-e2*np.sin(lat)*np.sin(lat))
    x=(v+h)*np.cos(lat)*np.cos(lon)
    y=(v+h)*np.cos(lat)*np.sin(lon)
    z=(v*(1-e2)+h)*np.sin(lat)

    return x,y,z



# (5) Gravitiy Models


# (5.1) Gravity Model ECEF - Frame
def Gravity_ECEF( r_eb_e ):
    #Gravitation_ECI - Calculates  acceleration due to gravity resolved about 
    #ECEF-frame
    #
    # Software for use with "Principles of GNSS, Inertial, and Multisensor
    # Integrated Navigation Systems," Second Edition.
    #
    # This function created 1/4/2012 by Paul Groves
    #
    # Inputs:
    #   r_eb_e  Cartesian position of body frame w.r.t. ECEF frame, resolved
    #           about ECEF-frame axes (m)
    # Outputs:
    #   g       Acceleration due to gravity (m/s^2)

    # Copyright 2012, Paul Groves
    # License: BSD; see license.txt for details

    # Modified for Python by: Felix Esser, E-Mail: s7feesse@uni-bonn.de


    #Parameters
    R_0 = 6378137 #WGS84 Equatorial radius in meters
    mu = 3.986004418E14 #WGS84 Earth gravitational constant (m^3 s^-2)
    J_2 = 1.082627E-3 #WGS84 Earth's second gravitational constant
    omega_ie = 7.292115E-5  #Earth rotation rate (rad/s)

    # Begins

    # Calculate distance from center of the Earth
    mag_r = np.linalg.norm(r_eb_e)

    g = np.array((0,0,0),dtype=float)

    # If the input position is 0,0,0, produce a dummy output
    if (mag_r == 0):
        return g
    else:
        z_scale = 5 * (r_eb_e[2] / mag_r)**2
        gamma = (-mu / mag_r**3) * (r_eb_e + 1.5 * J_2 * (R_0 / mag_r)**2 * np.array(((1 - z_scale) * r_eb_e[0], (1 - z_scale) * r_eb_e[1], (3 - z_scale) * r_eb_e[2])))
        # Add centripetal acceleration using (2.133)

        g[0] = gamma[0] + omega_ie**2 * r_eb_e[0]
        g[1] = gamma[1] + omega_ie**2 * r_eb_e[1]
        g[2] = gamma[2]
    return g


def Mat2TriUpright( M ):
    # matrix (must be symmetric! ) to upper triangle matrix

    # extract upper Triangle matrix
    return M[np.triu_indices( len(M[:,0]) )] 






