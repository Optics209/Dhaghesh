import functions.geodetic_tools as geod
import matplotlib.pyplot as plt
import numpy as np
import functions.geodetic_tools as geod

# Class for reading IMAR data from .mat file 
# @author: Felix Esser, s7feesse@uni-bonn.de
# @date: 14.10.20

# data structure:
# 0    # 1   # 2    # 3   # 4          # 5         # 6         # 7    # 8     # 9   #
# GPST # Lat # Long # alt # velo north # velo east # velo down # roll # pitch # yaw # rpy in [deg] !

class StrapdownData:
    def __init__(self):
        self.gpstime = np.empty( (0,1), dtype = float) 
        self.lat = np.empty( (0,1), dtype = float) 
        self.lon = np.empty( (0,1), dtype = float) 
        self.alt = np.empty( (0,1), dtype = float) 
        self.velo = np.empty( (0,3), dtype = float) 
        self.rpy = np.empty( (0,3), dtype = float)
        self.UTM = np.empty( (0,3), dtype = float)


    def read(self, filelocation, flag  ):

        if flag == "advanced":
            
            strap_data = np.loadtxt( filelocation, dtype=np.float)

            self.gpstime = strap_data[:,0]
            self.lat = strap_data[:,1]
            self.lon = strap_data[:,2]
            self.alt = strap_data[:,3]
            self.velo = strap_data[:,4:7]
            self.rpy = strap_data[:,7:10]

            # convert to UTM coordinates
            # UTM and height (own solution)
            N, E = geod.ell2utm( geod.deg2rad(self.lat), geod.deg2rad(self.lon), 32 * np.ones( self.lon.shape ) )
            self.UTM = np.c_[E,N,self.alt]
        else:
            strap_data = np.loadtxt( filelocation, dtype=np.float, delimiter = ",")

            self.gpstime = strap_data[:,0]
            self.lat = strap_data[:,1]
            self.lon = strap_data[:,2]
            self.alt = strap_data[:,3]
            self.velo = strap_data[:,4:7]
            self.rpy = strap_data[:,7:10]

            # convert to UTM coordinates
            # UTM and height (own solution)
            N, E = geod.ell2utm( geod.deg2rad(self.lat), geod.deg2rad(self.lon), 32 * np.ones( self.lon.shape ) )
            self.UTM = np.c_[E,N,self.alt]



