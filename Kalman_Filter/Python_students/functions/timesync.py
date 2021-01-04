import numpy as np



def timesync(gps_time, imu_time, thres):

    # offset between GPS and IMU
    idx = np.where( (gps_time - imu_time[0]) < 0.0 )

    # initialization index variables
    idx_imu = np.empty((1,0), dtype = int) 
    idx_gps = np.empty((1,0), dtype = int) 

    # gps counter
    gps_j = len(idx[0])

    # threshold
    thres = 0.0005

    for i in range(0, len(imu_time) ):
        
        # absolute time stamp difference
        delta_t = np.absolute(gps_time[gps_j] -  imu_time[i])

        # Intersection point
        if ( delta_t < thres ):
            # save indices
            idx_imu = np.append(idx_imu, i)
            idx_gps = np.append(idx_gps, gps_j)

            # update GPS counter
            gps_j +=1

        # stop if end of GPS is reached   
        if ( gps_j >= len(gps_time) ):
            break
        
    return idx_gps, idx_imu