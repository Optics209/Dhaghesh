import numpy as np
import scipy.io
from scipy import interpolate
from numpy.linalg import inv
import matplotlib.pyplot as plt


# from functions.geodetic_tools import * 
from functions.prediction import prediction
from functions.timesync import timesync
from functions.readIMARdata import *


## Mobile Sensing and Robotics - Exercise 2 - IGG Bonn, 12.10.20
# Author 1: M.Sc. Erik Heinz,
# Author 2: M.Sc. Tomislav Medic, 
# Author 3: B.Sc. Felix Esser (Python Code)
# Contact: medic@igg.uni-bonn.de, Tomislav Medic
#          s7feesse@uni-bonn.de, Felix Esser (Python)

# -----------------------------------------------------
# read data
imar_data = IMARdata()
imar_data.readIMAR("data/IMAR.mat", correctedIMUdata = True)

# -----------------------------------------------------
# Measurements

# Info IMU
f = 1000                  # Measuring frequency [Hz]
dt = 1/f                  # Measuring rate [s]

# Measurements
x = imar_data.gpsUTM[:,0]                 # East UTM [m]
y = imar_data.gpsUTM[:,1]                 # North UTM [m]
a = -imar_data.acceleration[:,0]          # acceleration x [m/s^2]
omega = imar_data.angularvelocity[:,2]    # angular velocity z [rad/s]

# -----------------------------------------------------
# Time Synchronisation
idx_gps, idx_imu = timesync( imar_data.gpstime, imar_data.imutime, 0.0005 )

# -----------------------------------------------------
# GPS Error Simulation

# (1) Simulation with normal distrubution
GPS_err_simu = False

if (GPS_err_simu == True):

    # index noize increase
    bd_idx = np.arange(start=700, stop=750, step=1, dtype=int)

    # parameter for normal distribution
    mu_e, sigma_e = 0, 1
    mu_n, sigma_n = 0, 1

    # values from normal distribution
    east_err = np.random.normal( mu_e, mu_e, len(bd_idx) )
    north_err = np.random.normal( mu_n, sigma_n, len(bd_idx) )

    # add noize to GPS observations
    x[bd_idx] += east_err
    y[bd_idx] += north_err

# (2) Multipath Simulation (systematic error)
GPS_multi_simu = False 

if (GPS_multi_simu == True):

    # length of GPS breakdown index
    bd_idx = np.arange(start=1200, stop=1250, step=1, dtype=int)

    # add multipath systematic to GPS observations 
    x[bd_idx] += 2
    y[bd_idx] += -1.3

# TODO (3) choose noise model for GPS

# --------------------------------------------------- #
# ########### Extended Kalman Filter ################ #
# --------------------------------------------------- #

# -----------------------------------------------------
# Statistics
print(np.std(omega[0:1000]))
# Measurement noise
std_gps = 0.05      # GPS [m]
std_a =  np.std(a[0:1000])       # Accelerations [m/s^2]
std_omega =  np.std(omega[0:1000])  # Angular rate [rad/s]

# System noise
wk_phi =  0.1   # Angular accelerations [rad/s]
wk_a =    1.5   # Linear jerk [m/s^3]

# # Covariance Matrix Measurements
Sll = np.array([[std_gps, 0, 0, 0],
                [0, std_gps, 0, 0],
                [0, 0, std_a,0, 0],
                [0, 0, 0, std_omega]])

# # Covariance Matrix System Noise
S_wkwk = np.array([[wk_phi, 0],
                   [0,   wk_a]])

# # Covariance Matgrix Initial States
# S_xkxk = 
# # -----------------------------------------------------


# -----------------------------------------------------
# Measurements

# define zero block
zero_block = np.nan * np.zeros( (2,len(imar_data.imutime) ) )

# Measurement matrix with IMU measurements
L = np.vstack(( np.transpose(imar_data.imutime) , zero_block, np.transpose(a), np.transpose(omega) ))
# add GPS at intersection time stamps
L[1:3,idx_imu] = np.vstack(( np.transpose( x[idx_gps] ), np.transpose( y[idx_gps] ))) 

# -----------------------------------------------------
# Initial states
# xk = 

# # -----------------------------------------------------
# # Variables to save estimated states

# # states
# xstate = np.zeros((len(imar_data.imutime),7))
# xstate[0,0] = imar_data.imutime[0]
# xstate[0,1:7] = xk
# print(xstate)
# # -------------------------------------------

# # Designmatrix (Jacobian)
# H = np.array([[1, 0, 0, 0, 0, 0],
#               [0, 1, 0, 0, 0, 0],
#               [0, 0, 0, 0, 0, 1],
#               [0, 0, 0, 1, 0, 0]])

# # -------------------------------------------
# # loop info
# percent = 0.05
# nbr = len(imar_data.imutime) # all = len(imar_data.imutime)
# print('filter started ... ')

# # -------------------------------------------
# # MAIN KF Loop

# for i in range(0, nbr ): 
#     # -------------------------------------------
#     # percent update 
#     if (i / nbr) > percent:
#         print('%1.0f' % (percent * 100), '% / 100%' )
#         percent += 0.05

#     # -------------------------------------------
#     # Prediction Step
#     x_bar, Sx_bar = prediction( xk, S_xkxk, S_wkwk, dt, i )

#     # -------------------------------------------
#     # Update Step (IMU only)

#     if ( np.isnan( L[1,i]) == True ):

#         # Kalman Gain Matrix
#         K = Sx_bar @ H[2:,:].T @ np.linalg.inv(H[2:,:] @ Sx_bar @ H[2:,:].T + Sll)

#         # Update state vector
#         x_dach = xk + K @ (L[-2:,i] - (H[-2:,:] @ x_bar))

#         # Covariance matrix states 
#         Sx_dach = (np.identity(Sx_bar.shape[0]) - (K @ H[-2:,:])) @ Sx_bar

#         # save current estimate
#         xstate[i,0] = imar_data.imutime[i]
#         xstate[i,1:7] = x_dach
        
#         # update states for next iteration
#         xk = x_dach
#         S_xkxk = Sx_dach
    
#     # -------------------------------------------
#     # Update (GPS + IMU)
#     else:

#         # Kalman Gain Matrix
#         K = Sx_bar @ H.T @ np.linalg.inv(H @ Sx_bar @ H.T + Sll)

#         # Update state vector
#         x_dach = xk + K @ (L[i] - (H @ x_bar))

#         # Covariance matrix states 
#         Sx_dach = (np.identity(Sx_bar.shape[0]) - (K @ H)) @ Sx_bar

#         # save estimation
#         xstate[i,0] = imar_data.imutime[i]
#         xstate[i,1:7] = x_dach

#         # update states for next iteration
#         xk = x_dach
#         S_xkxk = Sx_dach

# # end MAIN loop
# # -------------------------------------------

# print( 100, '%')
# print('... done')

# # --------------------------------------------------- #
# # ########### Plot EKF results ################ #
# # --------------------------------------------------- #

# # -------------------------------------------
# # reduce plot by indexing
# idx_plot = np.arange(start=0, stop=nbr, step=10, dtype=int)

# # -------------------------------------------
# # Trajectory plot
# plt.plot(x,y, '.b', markersize=12)
# plt.plot(xstate[idx_plot,1], xstate[idx_plot,2], '.r')
# plt.axis('equal')
# plt.title('EKF Trajectory and GPS Measurements', fontsize=14, fontweight='bold')
# plt.xlabel('UTM (East) [m]', fontsize=12, fontweight='bold')
# plt.ylabel('UTM (North) [m]', fontsize=12, fontweight='bold')
# plt.legend(['GPS Measurements', 'EKF Trajectory'])
# plt.grid(color='k', linestyle='-', linewidth=0.5)
# plt.show()

# # -------------------------------------------
# # Acceleration & Angular Velocity

# plt.subplot(211)
# plt.plot( imar_data.imutime , a, '.b' )
# plt.plot( xstate[idx_plot,0], xstate[idx_plot,6], '-r' )
# plt.ylabel(" Acceleration [m/s^2] ", fontsize=12, fontweight='bold')
# plt.legend(["Raw IMU Accelerations", "EKF Accelerations "])
# plt.title("IMU Accelerations vs. Filtered Accelerations (EKF) [X-ACC]", fontsize=14, fontweight='bold' )
# plt.grid(color='k', linestyle='-', linewidth=0.5)
    
# plt.subplot(212)
# plt.plot( imar_data.imutime , geod.rad2deg(omega), '.b' )
# plt.plot( xstate[idx_plot,0], geod.rad2deg(xstate[idx_plot,4]), '-r' )

# plt.ylabel(" Angular Velocity [deg/s] ", fontsize=12, fontweight='bold')
# plt.xlabel(" seconds of day [s] ", fontsize=12, fontweight='bold')
# plt.legend(["Raw IMU Angular Velocity", "EKF Angular Velocity "])
# plt.title("IMU Accelerations vs. Filtered Accelerations (EKF) [YAW] ", fontsize=14, fontweight='bold' )
# plt.grid(color='k', linestyle='-', linewidth=0.5)
# plt.show()

# # -------------------------------------------

# plt.show()