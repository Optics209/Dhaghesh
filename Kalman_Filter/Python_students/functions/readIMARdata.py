import scipy.io
import numpy as np
import functions.geodetic_tools as geod

# Class for reading IMAR data from .mat file 
# @author: Felix Esser, s7feesse@uni-bonn.de
# @date: 14.10.20


class IMARdata:
    def __init__(self):
        
        # GPS variables 
        self.gps = [] # IMAR GPS data (list structure)
        self.gpstime = np.empty( (0,1), dtype = float)
        self.gpsLLA = np.empty( (0,3), dtype = float)
        self.gpsUTM = np.empty( (0,3), dtype = float)

        # IMU variables
        self.imu = [] # IMAR IMU data (list structure)
        self.imutime = np.empty( (1,0), dtype = float)
        self.acceleration = np.empty( (0,3), dtype = float)
        self.angularvelocity = np.empty( (0,3), dtype = float)
        self.rpy_ned = np.empty( (0,3), dtype = float)

    # function to read data from IMAR.mat file (MATLAB)
    def readIMAR(self, filelocation, correctedIMUdata ):

        data = scipy.io.loadmat( filelocation )

        # GPS Observations and IMU Measurements
        self.gps = data['gps']
        gps_time = self.gps['gps_time'][0,0]
        self.gpstime = gps_time[:,0]
        self.gpsLLA = self.gps['gps_xyz'][0,0] # LLA !

        # calclulate UTM Coordinates
        N,E = geod.ell2utm( self.gpsLLA[:,1], self.gpsLLA[:,0], 32 * np.ones(len(self.gpsLLA[:,0])) )

        self.gpsUTM = np.c_[E,N,self.gpsLLA[:,2]]

        # IMU raw data
        self.imu = data['imu']
        imu_time = self.imu['imu_time'][0,0]
        self.imutime = imu_time[:,0]
        self.rpy_ned = self.imu['rpy_ned'][0,0]

        if ( correctedIMUdata == True ):
            self.acceleration   = self.imu['acc_corr'][0,0]
            self.angularvelocity = self.imu['omg_corr'][0,0]
        else:
            self.acceleration   = self.imu['acc_ib_b'][0,0]
            self.angularvelocity = self.imu['omg_ib_b'][0,0]








	