import numpy as np
import matplotlib.pyplot as plt

import scipy.io
from scipy import interpolate

from functions.readIMARdata import *
import functions.geodetic_tools as geod

#from functions.strapdown_functions import *


## Determine where your m-file's folder is.
#folder = fileparts(which(mfilename)); 
#Add that folder plus all subfolders to the path.
#addpath(genpath(folder));
# -------------------------------------------------------------------------
#                     STRAPDOWN NAVIGATION 

# for: Master of Science in Geodetic Engineering
#      I. semester - Sensors & State Estimation
#      Geodetic track - Exercise 1: Strapdown navigation

# -------------------------------------------------------------------------
# Description: 
# * Stand-alone solution
# * Made for processing IMAR data file from IGG kinematic MSS system
# * Navigation equations in local navigation frame (NED)
# * Algorithm reference: Groves 2013. "Principles of GNSS, Inertial, and 
#   Multisensor Integrated Navigation Systems," Second Edition, Chapter 5.
# * All symbols/variables are defined corresponding to the ref. literature
# * Summary of all relevant functions can be found in:
#   "MGE-MSR-01 Exercise 1 - Supplementary Document"
# -------------------------------------------------------------------------
# Author1: Tomislav Medic, 14.11.2017. IGG Bonn ( MATLAB Implementation )
# Author2: Felix Esser, 09.10,2020. IGG Bonn ( PYTHON Implementation )



# #########################################################################
# CONTROL VARIABLES
# ##########################################################################

# Simplification:
# (simplification means that Earths rotation rate, transportation rate, and
# gravity influence in North direction are disregarded)
simplified = 0 # 0 - no, 1 - yes  

# Empirical gravity model:
# Instead of mathematical gravity model uses accelerometer readings at
# still-state to estimate the acceleration due to gravity 
#(1st 1 second of the trajectory - 1000 observations)
empirical_gravity = 0
obs_nr = 40000 # length emp
# Trajectory analysis (visual & numerical analysis)
t_analysis = 1 # 1 - yes, 0 - no
# Saving navigation profile (out_profile)
save_solution = 0 # 1 - yes, 0 - no


# ##########################################################################
#  CONSTANTS ,HELP VARIABLES, FUNCTIONS
# ##########################################################################

# Earth's ellipsoid radius at equator - a [m]
a = 6378137; #GRS80 & WGS84
# Earth's ellipsoid radius at poles - b [m]
# b = 6356752.31425; % WGS84
b = 6356752.31414 # GRS80
# Earth's ellipsoid eccentricity - e 
# e = 0.0818191908425; % WGS84
e = 0.0818191910428 # GRS80
# Earth's rotational speed [rad s^-1]
w_ie = 7.2921150 * 10**-5 # 7.292115*10**-5
w_ie_vec = np.array((0, 0, w_ie)) # in ECEF (GRS80) coordinates

# Geocentric gravitational constant - mi [m^3 s^-2]
mi = 3.986004418 * 10**14

# Earth's ellipsoid flattening - f
# f = 1 / 298.257223563; % WGS84
f = 1 / 298.257222101 #GRS80

# ##########################################################################
# 0. Loading Data
# ##########################################################################

# adapted for: IMAR data file from IGG kinematic MSS (Erik's Time Machine)

# load path and filename
#data = scipy.io.loadmat('../data/IMAR.mat')

imar_data = IMARdata()
#imar_data.readIMAR('../data/IMAR0007_cut01_matlab_sparse.mat', correctedIMUdata = False, UTMZone=32, fromMatlab=True)
imar_data.readIMAR('../data/IMAR.mat', correctedIMUdata = False )

# --------------------------------------------------------------
# IMU DATA

# Accelerometer measurements (raw)
# matrix containing all measured accelerations f_ib_b = [f_x,f_y,f_z] (nx3) [m s^-2]
allmsm_f_ib_b = imar_data.acceleration

# Gyro measurements (raw)
# matrix containing all measured angular rates in b-frame w_ib_b=[w_x,w_y,w_z] (nx3)[rad s^-1]
allmsm_w_ib_b = imar_data.angularvelocity

# Roll Pitch Yaw in NED Frame (Navigation frame)
rpy_ned = imar_data.rpy_ned

# Time stamps 
all_time = imar_data.imutime

# Number of measurement epochs (n)
no_epochs = len(all_time)

# --------------------------------------------------------------
# GPS DATA
gps_xyz = imar_data.gpsLLA
gps_time = imar_data.gpstime

# --------------------------------------------------------------
# Navigation solution/profiles matrix
# Data format of output matrix [no_epochs x 10]:
# [ Time | Position (1x3) | Velocity (1x3) | Attitude (1x3) ]
#   [s]  |     [m,m,m]    |    [m s^-1]    |      [deg]

nav_solution = np.zeros( (no_epochs, 10) ) 

# ##########################################################################
# 1. Defining Initial Values
# ##########################################################################

# 1.0. Initial Time (t [s])

# (second of the day synchronized for IMU and GNSS)
time_old = all_time[0]

# 1.1. Initial Position (Lat, Long, h [rad,rad,m])
# Taken from GNSS observations of the MSS (Mulitsensor System)
if ( gps_xyz[0,0] < (2*np.pi) ):
    Lat_b_old = gps_xyz[0,1]
    Long_b_old = gps_xyz[0,0]
    h_b_old = gps_xyz[0,2]

    # ECEF Initial Position:
    xyz0 = geod.ell2xyz( Lat_b_old, Long_b_old, h_b_old )
    LLA0 = geod.vector3( Lat_b_old, Long_b_old, h_b_old )

else:
    xx = gps_xyz[0,0]
    yy = gps_xyz[0,1]
    zz = gps_xyz[0,2]
    Lat_b_old, Long_b_old, h_b_old = geod.ecef2lla( xx, yy, zz, ellipsoid= "GRS80", version=1 )

    # ECEF Initial Position
    xyz0 = geod.vector3(xx,yy,zz)
    LLA0 = geod.vector3( Lat_b_old, Long_b_old, h_b_old )


# 1.2 Initial attitude/coordinate transformation matrix CTM 
# Attitude - Euler angles - roll pitch yaw/heading 
# NED w.r.t. ECEF [rad] [rad] [rad]

eul_n_b_old = rpy_ned[0,:] # start Euler angles from IMAR solution

# ---------------------------------------
# TODO (1) INITIAL VALUES

# - set up rotation matrix here!

C_b_n_old = 

# 1.3. Initial velocity (North, East, Down [m s^-1])
v_eb_n_old = 

# ---------------------------------------

# 1.4 Storing initial navigation solution
nav_solution[0,0] = time_old
nav_solution[0,1] = geod.rad2deg(Lat_b_old)
nav_solution[0,2] = geod.rad2deg(Long_b_old)
nav_solution[0,3] = h_b_old
nav_solution[0,4:7] = v_eb_n_old
nav_solution[0,7:10] = geod.rad2deg(eul_n_b_old)

# Calculate empirical gravity in B-Frame
# TODO Compute empirical Gravity from IMU measurements (eq. 16)
g_emp = 

# percent update variable
percent = 0.05

# ----------------------------------------------------------------------------

# 1.5. Start of the main loop
for epoch in range(1, no_epochs):  # for testing: 50000 # all no_epochs
    
    # percent update
    if (epoch / no_epochs) > percent:
        print('%1.0f' % (percent * 100), '% / 100%' )
        percent += 0.05
        
    # 1.6. Time interval of epoch i, where i = 1,2..no_epochs - tau [s]
    time_new = all_time[epoch]
    tau_i = time_new - time_old
    #tau_i = tau_i[0] ##array to float
    
    # 1.7. Measurements of epoch i 
    # accelerometer measurements - specific force/accelaration
    f_ib_b = allmsm_f_ib_b[epoch-1,:]
    # gyro measurements - angular rates
    w_ib_b = allmsm_w_ib_b[epoch-1,:]

## ------------------------------------------------------------------- ##
## STRAPDOWN N FRAME 
## ------------------------------------------------------------------- ##

## 2. Attitude update (Coordinate transformation matrix - CTM)
# ------------------------------------------------------------------------

    # Skew symmetric matrix of gyro measurements
    # ---------------------------------------
    # TODO (2) ATTITUDE, equation 3
    # - set up skew symmetric matrix of measured angular velocity

    OM_ib_b = 
            
    # Matrix for Earth's rotation rate
    # ---------------------------------------
    # TODO (2) ATTITUDE, equation 3
    # - set up matrix that represents earth rotattion rate

    OM_ie_n = 

    # Simple CTM update (simplifications)
    # ...................................
 
    if (simplified == 1):

        # CTM update, simplified
        # TODO (2) ATTITUDE
        # - perform Attitude Update with simplified equation

        C_b_n_new = 
       
    # Regular CTM update (no simplifications)
    # .......................................
    else:
        # Matrix for Transport rate

        # Curvature radius R_e in the East direction
        # TODO 
        # - Re (eq. 8)
        R_e = 
     
        # Curvature radius R_n in the North direction
        # TODO 
        # - Rn (eq. 9)
        R_n = 

        # Vector of angular rates of NED w.r.t. ECEF
        # TODO 
        # - compute velocity of body frame (eq. 7)
        w_en_n = 

        # Skew symmetric matrix for Transport rate
        # TODO 
        # - set up skew symmetric Matrix of Transport rate
        OM_en_n =

        # CTM update
        # ..................................................................
        # TODO (2) ATTITUDE
        # - perform Attitude Update, not simplified
        C_b_n_new = 

    # 3. Specific force frame transofrmation
    # ---------------------------------------------------------------------

    # TODO (3) SPECIFIC FORCE TRANSFORMATION
    # - transform specific force to Navigation Frame (eq. 10)
    f_ib_n = 

    ## 4. Velocity update
    # --------------------------------------------------------------------- 
    
    # use same variable for empirical and model gravity
    g_b_n = geod.vector3( 0, 0, 0 )

    # model gravity
    # TODO (4) VELOCITY UPDATE
    # - compute model gravity from model (eq. 15)
    # - D --> Down component
    
    if (empirical_gravity == 0):
        # TODO (4) VELOCITY UPDATE
        # - compute gravity from model
        g_0 = 
        g_b_n[2] = 
    else:
        # Gravity calculated based on accelerometer measurements
        # TODO (4) VELOCITY UPDATE
        # - use empirical gravity, computed from measurements
        # - hint: compute outside loop
        g_b_n[2] = g_emp
    
    # Simplified: without g_b_n in north direction & transport rate
    if (simplified == 1):
        
        # TODO (4) VELOCITY UPDATE
        # - implement simplified version of velocity update (eq. 11)
        v_eb_n_new =

    else:
        # TODO (4) VELOCITY UPDATE
        # - compute North Component of Gravity Vector (eq. 14)
        g_b_n[0] = 

        # TODO (4) VELOCITY UPDATE
        # - perform not simplified velocity update (eq. 11)
        v_eb_n_new = 
    
    # (5) Position Update
    # ---------------------------------------------------------------------
    
    # Updating ellipsoid height
    # TODO (5) Position UPDATE
    # - perform height update (eq. 17)
    h_b_new = 
    
    if( simplified == 1):
       # TODO (5) Position UPDATE
       # - implement simplified version
       # Curvature radius R_e in the East direction
       R_e = 
       
       # Curvature radius R_n in the North direction

       R_n = 
    
    # Updating lattitude
    # TODO (5) Position UPDATE
    # - compute latitude update (eq. 18)
    Lat_b_new = 
    
    # Calculating new curvature radius in the East direction (R_e_new)
    # TODO (5) Position UPDATE
    # - compute new curvature in east direction --> needed for longitude update!
    R_e_new = 
    
    # Updating longitude
    # TODO (5) Position UPDATE
    # - compute longitude update (eq. 19)
    Long_b_new = 

    # 5. Storing values in navigation solution matrix
    # ---------------------------------------------------------------------    
    nav_solution[epoch,0] = time_new
    nav_solution[epoch,1] = geod.rad2deg(Lat_b_new)
    nav_solution[epoch,2] = geod.rad2deg(Long_b_new)
    nav_solution[epoch,3] = h_b_new
    nav_solution[epoch,4:7] = v_eb_n_new

    # Rotation Matrix to Euler Angles conversion
    rpy = geod.Rotmat2Euler(C_b_n_new)
    
    # save Orientation to Matrix
    nav_solution[epoch,7:10] = geod.rad2deg( rpy )
    
    # 6. Updating states
    # ---------------------------------------------------------------------    
    Lat_b_old = Lat_b_new
    Long_b_old = Long_b_new
    h_b_old = h_b_new
    v_eb_n_old = v_eb_n_new
    C_b_n_old = C_b_n_new
    time_old = time_new

print('100 % ')

# Data Structure
# [ Time | Position (1x3) | Velocity (1x3) | Attitude (1x3) ]
#   [s]  |     [m,m,m]    |    [m s^-1]    |      [deg]

# save and load Strapdown Estimates
save_solution == True

if (save_solution == True):
    np.savetxt('nav_solution_python.txt', nav_solution, delimiter=',')
    print('saved nav solution file')
    # nav_solution = np.loadtxt( 'nav_solution_python.txt', dtype=np.float, delimiter=',' )

# ######################################################################### 
                    # ANALYSIS OF THE RESULTS 
# ######################################################################### 

# TODO Open ANALYSIS TOOL