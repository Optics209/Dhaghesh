# Analysis Tool Strapdown

import numpy as np
import matplotlib.pyplot as plt

import scipy.io
from scipy import interpolate

import functions.geodetic_tools as geod
from functions.readIMARdata import * 
from functions.readStrapdownData import *


# (1) Load reference data (Grooves fancy strapdown)
str_adv = 'analysisTool/data/Groves_accurate_solution.txt'
SD_advanced = StrapdownData()
SD_advanced.read( str_adv, "advanced")

# (2) Load Own Solution

str_simple = 'analysisTool/data/nav_solution_python.txt'
SD_simple = StrapdownData()
SD_simple.read( str_simple , "simple")

# (3) Impirical Gravity One

str_emprical = 'analysisTool/data/nav_solution_python_emp_gravity_one.txt'
SD_emprical = StrapdownData()
SD_emprical.read( str_emprical , "empirical gravity")

# (4) load reference data
imar = IMARdata()
imar.readIMAR('analysisTool/data/IMAR_2018.mat', correctedIMUdata = False )


# INTERPOLATE DATA

# time delta of interpolation
dt = 0.1

# start and end time
start_time = np.amax( geod.vector3(imar.gpstime[0], SD_simple.gpstime[0], SD_advanced.gpstime[0]) )
end_time = np.amin( geod.vector3(imar.gpstime[-1], SD_simple.gpstime[-1], SD_advanced.gpstime[-1]) )

# New GPS time vector
time_vector = np.arange(start=start_time, stop=end_time, step=dt, dtype=float)
time_vector_red = time_vector - time_vector[0]

# (1) INTERPOLATION UTM and Heights
# easting
f_e_gps = interpolate.interp1d( imar.gpstime, imar.gpsUTM[:,0], kind='linear')
E_GPS = f_e_gps(time_vector)

f_e_imu = interpolate.interp1d( SD_simple.gpstime, SD_simple.UTM[:,0], kind='linear')
E_sim = f_e_imu(time_vector)

f_e_imu = interpolate.interp1d(  SD_advanced.gpstime, SD_advanced.UTM[:,0], kind='linear')
E_adv = f_e_imu(time_vector)

# Northing
f_n_gps = interpolate.interp1d( imar.gpstime, imar.gpsUTM[:,1], kind='linear')
N_GPS = f_n_gps(time_vector)

f_n_imu = interpolate.interp1d( SD_simple.gpstime, SD_simple.UTM[:,1] ,kind='linear')
N_sim = f_n_imu(time_vector)

f_n_imu = interpolate.interp1d(  SD_advanced.gpstime, SD_advanced.UTM[:,1] ,kind='linear')
N_adv = f_n_imu(time_vector)

# heights
f_n_gps = interpolate.interp1d( imar.gpstime, imar.gpsUTM[:,2], kind='linear')
H_GPS = f_n_gps(time_vector)

f_n_imu = interpolate.interp1d( SD_simple.gpstime, SD_simple.UTM[:,2] ,kind='linear')
H_sim = f_n_imu(time_vector)

f_n_imu = interpolate.interp1d( SD_advanced.gpstime, SD_advanced.UTM[:,2] ,kind='linear')
H_adv = f_n_imu(time_vector)

f_n_imu = interpolate.interp1d( SD_emprical.gpstime, SD_emprical.UTM[:,2] ,kind='linear')
H_emp = f_n_imu(time_vector)

# (2) roll pitch yaw


# start and end time
start_time = np.amax( geod.vector3(imar.imutime[0], SD_simple.gpstime[0], SD_advanced.gpstime[0]) )
end_time = np.amin( geod.vector3(imar.imutime[-1], SD_simple.gpstime[-1], SD_advanced.gpstime[-1]) )

# New GPS time vector
time_vector_rpy = np.arange(start=start_time, stop=end_time, step=dt, dtype=float)
time_vector_red_rpy = time_vector - time_vector[0]


# IMAR
f_n_roll = interpolate.interp1d( imar.imutime, imar.rpy_ned[:,0], kind='linear')
f_n_pitch = interpolate.interp1d( imar.imutime, imar.rpy_ned[:,1], kind='linear')
f_n_yaw = interpolate.interp1d( imar.imutime, imar.rpy_ned[:,2], kind='linear')

roll_imar = f_n_roll(time_vector)
pitch_imar = f_n_pitch(time_vector)
yaw_imar = f_n_yaw(time_vector)

# SIMPLE
f_n_roll = interpolate.interp1d( SD_simple.gpstime, SD_simple.rpy[:,0], kind='linear')
f_n_pitch = interpolate.interp1d( SD_simple.gpstime, SD_simple.rpy[:,1], kind='linear')
f_n_yaw = interpolate.interp1d( SD_simple.gpstime, SD_simple.rpy[:,2], kind='linear')

roll_sim = f_n_roll(time_vector)
pitch_sim = f_n_pitch(time_vector)
yaw_sim = f_n_yaw(time_vector)

# ADVANCED
f_n_roll = interpolate.interp1d( SD_advanced.gpstime, SD_advanced.rpy[:,0], kind='linear')
f_n_pitch = interpolate.interp1d( SD_advanced.gpstime, SD_advanced.rpy[:,1], kind='linear')
f_n_yaw = interpolate.interp1d( SD_advanced.gpstime, SD_advanced.rpy[:,2], kind='linear')

roll_adv = f_n_roll(time_vector)
pitch_adv = f_n_pitch(time_vector)
yaw_adv = f_n_yaw(time_vector)


# PLOT

plotvar = True

if (plotvar == True):

    # (1) PLOT Trajectory
    plt.plot(E_sim, N_sim, '.r')
    plt.plot(E_adv, N_adv, ".g")
    plt.plot(E_GPS, N_GPS, '.k')

    plt.title('UTM Trajectory (GPS / Strapdown) interpolated', fontsize=14, fontweight='bold' )
    plt.axis('equal')
    plt.xlabel(" Easting [m] ", fontsize=12, fontweight='bold')
    plt.ylabel(" Nothing [m] ", fontsize=12, fontweight='bold')
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.legend(['Strapdown simple','Strapdown advanced', "GPS observations"])
    plt.show()

    # (2) PLOT Coordinate Differences
    plt.subplot(311)
    plt.plot(time_vector_red, E_sim - E_GPS, '.r')
    plt.plot(time_vector_red, E_adv - E_GPS, '.g')
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.title('UTM Coordinates - Differences to GPS observations', fontsize=14, fontweight='bold' )
    plt.legend(['simple solution - GPS', "advanced solution - GPS"], fontsize=10)
    plt.ylabel("$\Delta$ UTM east [m]", fontsize=12, fontweight='bold')

    plt.subplot(312)
    plt.plot(time_vector_red, N_sim - N_GPS, '.r')
    plt.plot(time_vector_red, N_adv - N_GPS, '.g')
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.legend(['simple solution - GPS', "advanced solution - GPS"], fontsize=10)
    plt.ylabel("$\Delta$ UTM north [m]", fontsize=12, fontweight='bold')

    plt.subplot(313)
    plt.plot(time_vector_red, H_sim - H_GPS, '.r')
    plt.plot(time_vector_red, H_adv - H_GPS, '.g')
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.legend(['simple solution - GPS', "advanced solution - GPS"], fontsize=10)
    plt.xlabel("time [s]", fontsize=12, fontweight='bold')
    plt.ylabel("$\Delta$ height [m]", fontsize=12, fontweight='bold')
    plt.show()

    # (3) PLOT RPY

    # RPY
    # IMAR
    plt.plot(time_vector_red_rpy, geod.rad2deg(roll_imar), '.r')
    plt.plot(time_vector_red_rpy, geod.rad2deg(pitch_imar), '.g')
    plt.plot(time_vector_red_rpy, geod.rad2deg(yaw_imar), '.b')

    plt.title('Roll Pitch Yaw IMAR', fontsize=14, fontweight='bold' )
    plt.xlabel(" time [s] ", fontsize=12, fontweight='bold')
    plt.ylabel(" [°] ", fontsize=12, fontweight='bold')
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.legend(['roll angle','pitch angle', "yaw angle"])
    plt.show()

    # Simple
    plt.plot(time_vector_red_rpy, roll_sim, '.r')
    plt.plot(time_vector_red_rpy, pitch_sim, '.g')
    plt.plot(time_vector_red_rpy, yaw_sim, '.b')

    plt.title('Roll Pitch Yaw Simple Strapdown', fontsize=14, fontweight='bold' )
    plt.xlabel(" time [s] ", fontsize=12, fontweight='bold')
    plt.ylabel(" [°] ", fontsize=12, fontweight='bold')
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.legend(['roll angle','pitch angle', "yaw angle"])
    plt.show()

    # advanced
    plt.plot(time_vector_red_rpy, roll_adv, '.r')
    plt.plot(time_vector_red_rpy, pitch_adv, '.g')
    plt.plot(time_vector_red_rpy, yaw_adv, '.b')

    plt.title('Roll Pitch Yaw Advanced Strapdown', fontsize=14, fontweight='bold' )
    plt.xlabel(" time [s] ", fontsize=12, fontweight='bold')
    plt.ylabel(" [°] ", fontsize=12, fontweight='bold')
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.legend(['roll angle','pitch angle', "yaw angle"])
    plt.show()

    #emprical

    plt.plot(time_vector_red_rpy, H_sim - H_GPS, '.r')
    plt.plot(time_vector_red_rpy, H_emp - H_GPS, '.b')

    plt.title('Heights Differences', fontsize=14, fontweight='bold' )
    plt.xlabel(" time [s] ", fontsize=12, fontweight='bold')
    plt.ylabel(" meters[m] ", fontsize=12, fontweight='bold')
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.legend(['Simple_Height','Empirical_Height'])
    plt.show()

    # RPY differences
    delta_roll_sim = roll_sim - geod.rad2deg(roll_imar)
    delta_pitch_sim = pitch_sim - geod.rad2deg(pitch_imar)
    delta_yaw_sim = yaw_sim - geod.rad2deg(yaw_imar)

    # remove outlier
    idx = np.where(delta_roll_sim > 5)
    if (idx[0].size > 0):
        delta_roll_sim = np.delete( delta_roll_sim, idx[0][0] )
        time_mod_roll_sim = np.delete( time_vector_red_rpy, idx[0][0] )
    else:
        time_mod_roll_sim = time_vector_red_rpy

    idx = np.where(delta_pitch_sim > 5)
    if (idx[0].size > 0):
        delta_pitch_sim = np.delete( delta_pitch_sim, idx[0][0] )
        time_mod_pitch_sim = np.delete( time_vector_red_rpy, idx[0][0] )
    else:
        time_mod_pitch_sim = time_vector_red_rpy

    idx = np.where(delta_yaw_sim > 5)
    if (idx[0].size > 0):
        delta_yaw_sim = np.delete( delta_yaw_sim, idx[0][0] )
        time_mod_yaw_sim = np.delete( time_vector_red_rpy, idx[0][0] )
    else:
        time_mod_yaw_sim = time_vector_red_rpy


    delta_roll_adv = roll_adv - geod.rad2deg(roll_imar)
    delta_pitch_adv = pitch_adv - geod.rad2deg(pitch_imar)
    delta_yaw_adv = yaw_adv - geod.rad2deg(yaw_imar)

    # remove outlier
    idx = np.where(delta_roll_adv > 5)
    if (idx[0].size > 0):
        delta_roll_adv = np.delete( delta_roll_adv, idx[0][0] )
        time_mod_roll_adv = np.delete( time_vector_red_rpy, idx[0][0] )
    else:
        time_mod_roll_adv = time_vector_red_rpy

    idx = np.where(delta_pitch_adv > 5)
    if (idx[0].size > 0):
        delta_pitch_adv = np.delete( delta_pitch_adv, idx[0][0] )
        time_mod_pitch_adv = np.delete( time_vector_red_rpy, idx[0][0] )
    else:
        time_mod_pitch_adv = time_vector_red_rpy

    idx = np.where(delta_yaw_adv > 5)
    if (idx[0].size > 0):
        delta_yaw_adv = np.delete( delta_yaw_adv, idx[0][0] )
        time_mod_yaw_adv = np.delete( time_vector_red_rpy, idx[0][0] )
    else:
        time_mod_yaw_adv = time_vector_red_rpy

    plt.subplot(211)
    plt.plot(time_mod_roll_sim, delta_roll_sim, '.r')
    plt.plot(time_mod_pitch_sim, delta_pitch_sim, '.g')
    plt.plot(time_mod_yaw_sim, delta_yaw_sim, '.b')
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.title('RPY Simple Strapdown - Differences to IMAR Solution', fontsize=14, fontweight='bold' )
    plt.legend(['$\Delta$ roll', "$\Delta$ pitch", "$\Delta$ yaw"], fontsize=10)
    plt.ylabel("$\Delta$ [°]", fontsize=12, fontweight='bold')

    plt.subplot(212)
    plt.plot(time_mod_roll_adv, delta_roll_adv, '.r')
    plt.plot(time_mod_pitch_adv, delta_pitch_adv, '.g')
    plt.plot(time_mod_yaw_adv, delta_yaw_adv, '.b')
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.title('RPY Advanced Strapdown - Differences to IMAR Solution', fontsize=14, fontweight='bold' )
    plt.legend(['$\Delta$ roll', "$\Delta$ pitch", "$\Delta$ yaw"], fontsize=10)
    plt.ylabel("$\Delta$ [°]", fontsize=12, fontweight='bold')
    plt.xlabel(" time [s] ", fontsize=12, fontweight='bold')
    plt.show()

    # VELOCITY PLOT
    # start and end time
    start_time = np.amax( geod.vector2( SD_simple.gpstime[0], SD_advanced.gpstime[0]) )
    end_time = np.amin( geod.vector2( SD_simple.gpstime[-1], SD_advanced.gpstime[-1]) )

    # New GPS time vector
    time_vector = np.arange( start=start_time, stop=end_time, step=dt, dtype=float )
    time_vector_red = time_vector - time_vector[0]

    # SIMPLE
    f_n_vx = interpolate.interp1d( SD_simple.gpstime, SD_simple.velo[:,0], kind='linear')
    f_n_vy = interpolate.interp1d( SD_simple.gpstime, SD_simple.velo[:,1], kind='linear')
    f_n_vz = interpolate.interp1d( SD_simple.gpstime, SD_simple.velo[:,2], kind='linear')

    vx_sim = f_n_vx(time_vector)
    vy_sim = f_n_vy(time_vector)
    vz_sim = f_n_vz(time_vector)

    # ADVANCED
    f_n_vx = interpolate.interp1d( SD_advanced.gpstime, SD_advanced.velo[:,0], kind='linear')
    f_n_vy = interpolate.interp1d( SD_advanced.gpstime, SD_advanced.velo[:,1], kind='linear')
    f_n_vz = interpolate.interp1d( SD_advanced.gpstime, SD_advanced.velo[:,2], kind='linear')

    vx_adv = f_n_vx(time_vector)
    vy_adv = f_n_vy(time_vector)
    vz_adv = f_n_vz(time_vector)


    # PLOT VELOCITY
    
    # simple
    plt.plot(time_vector_red, vx_sim, '.r')
    plt.plot(time_vector_red, vy_sim, '.g')
    plt.plot(time_vector_red, vz_sim, '.b')

    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.title('Velocity Simple Strapdown', fontsize=14, fontweight='bold' )
    plt.legend(['$v_x$', "$v_y$", "$v_z$"], fontsize=10)
    plt.ylabel("$m_s^2$", fontsize=12, fontweight='bold')
    plt.xlabel("time [s]", fontsize=12, fontweight='bold')
    plt.show()

    # advanced
    plt.plot(time_vector_red, vx_adv, '.r')
    plt.plot(time_vector_red, vy_adv, '.g')
    plt.plot(time_vector_red, vz_adv, '.b')

    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.title('Velocity Advanced Strapdown', fontsize=14, fontweight='bold' )
    plt.legend(['$v_x$', "$v_y$", "$v_z$"], fontsize=10)
    plt.ylabel("$m_s^2$", fontsize=12, fontweight='bold')
    plt.xlabel("time [s]", fontsize=12, fontweight='bold')
    plt.show()

    # plot differences
    plt.plot(time_vector_red, vx_adv - vx_sim, '.r')
    plt.plot(time_vector_red, vy_adv - vy_sim, '.g')
    plt.plot(time_vector_red, vz_adv - vz_sim, '.b')

    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.title(' Velocity Differences (Simple - Advanced) Strapdown', fontsize=14, fontweight='bold' )
    plt.legend(['$v_x$', "$v_y$", "$v_z$"], fontsize=10)
    plt.ylabel("$m_s^2$", fontsize=12, fontweight='bold')
    plt.xlabel("time [s]", fontsize=12, fontweight='bold')
    plt.show()