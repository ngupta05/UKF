# UKF
Unscented Kalman Filter

This repo implements unscented Kalman filter for Udacity's Self Driving Car Nanodegree project.

Brief explanation of what each source file does:

src/main.cpp: Reads in input file line by line, creates Lidar/Radar update and call UKF's processLidar/Radar functions.
Saves the updated state after each step, to compute RMSE later
UKF.h/cpp: Implements Unscented Kalman Filter logic and also does fusion between Radar and Lidar updates
Tools.h/cpp: Implements RMSE calculator
Update.h: Defines datastructures for Lidar and Radar measurements

Eigen: Matrix library taken from https://github.com/udacity/CarND-Unscented-Kalman-Filter-Project/tree/master/src/Eigen

Build and run instructions:

Download the repo and 'cd' to main repo folder
mkdir build && cd build
cmake .. && make
Run: ./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt output.txt
