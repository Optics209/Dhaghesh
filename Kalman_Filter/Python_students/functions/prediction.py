# prediction function for Kalman Filter (2D)

import numpy as np
import matplotlib.pyplot as plt


def vector6(x1,x2,x3,x4,x5,x6):
    return np.array((x1,x2,x3,x4,x5,x6), dtype=float)  


def prediction( xk, S_xkxk, S_wkwk, dt, i ):
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Function to predict the state vector of the actual epoch k+1 in a Kalman
    # filter. The motion model from Example 2 in the lectures is used.
    #
    # Input:     xk     ... estimated state vector (epoch k)
    #            S_xkxk ... estimated covariance matrix (epoch k)
    #            S_wkwk ... system noise
    #            dt     ... sampling rate
    #
    # Output:    x_bar  ... pedicted state vector (epoch k+1)
    #            Sx_bar ... predicted covariance matrix (epoch k+1)
    #
    # Author:    M.Sc. Erik Heinz, M.Sc. Tomislav Medic
    # Contact:   heinz@igg.uni-bonn.de, medic@igg.uni-bonn.de
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #%
    # Predicted state vector
	
    x_bar = np.array([[xk[0,0] + np.cos(xk[2,0]) * (xk[4,0] * dt)],
                      [xk[1,0] + np.sin(xk[2,0]) * (xk[4,0] * dt)],
                      [xk[2,0] + xk[2,0] * dt],
                      [xk[3,0] + dt],
                      [xk[4,0] + xk[5,0] * dt],
                      [xk[5,0] + dt]])
   
    # Transition matrix
    Tk = np.array([[1, 0, -np.sin(xk[2,0])* (xk[4,0] * dt), 0, np.cos(xk[2,0]) * dt, 0],
                   [0, 1,  np.cos(xk[2,0])* (xk[4,0] * dt), 0, np.sin(xk[2,0]) * dt, 0],
                   [0, 0, 1, dt, 0, 0],
                   [0, 0, 0, 1,  0, 0],
                   [0, 0, 0, 0, 1, dt],
                   [0, 0, 0, 0, 0, 1]])
              
    # System noise matrix
    Sk = np.array([[0,0],
                   [0,0],
                   [0,0],
                   [dt,0],
                   [0,0],
                   [0,dt]])

    # Pedicted covariance matrix
    # Sx_bar = Tk * S_xkxk * Tk' + Sk * S_wkwk * Sk'
    Sx_bar = (Tk @ S_xkxk @ Tk.T) + (Sk @ S_wkwk @ Sk.T) 

    return x_bar, Sx_bar