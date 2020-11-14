# **Ego-vehicle Mathematical Modelling and State Estimation**

## Project goal
Implementation of a Race Car states estimation algorithm of Extended Kalman Filter (EKF) in MATLAB/Simulink environment. Stages of building a mathematical model of the vehicle equipped with state estimator are considered to achieve sensor fusion goal of compensating sensor errors. This project is a part of my Master's Thesis work.

## Project steps
1. Vehicle Mathematical Modelling: Bicycle model (roll and pitch rates neglected) to obtain ground truth data of the vehicle.

    <img src="./imgs/1_Bicycle model diagram.PNG" alt="Bicycle model">
    
    _Source: A. Katriniok et al., IEEE European Control Conference (ECC), Malaysia, 2013, - P. 975._
    
    **Simulink model**
    <img src="./imgs/2_Simulink single-track model.PNG" alt="Simulink bicycle model" width="800" height="550"> 

2. Sensors Modelling: Imitation of GPS and IMU sensors noisy measurement by adding Gaussian noise to the ground truth data.
    **Noisy GPS measurement of position coordinate x**
    <img src="./imgs/3_Coordinate x plot GPS sensor.png" alt="GPS x measurement" width="800" height="430"> 
    
3. Process Modelling: Representing vehicle motion using two non-linear process models of Constant Turn Rate and Velocity **CTRV** and  Constant Turn Rate and Acceleration **CTRA**
4. State Estimation: GPS and IMU noisy measurement are fused using Extended Kalman Filter (EKF) algorithm.

<img src="./imgs/4_Full System model.PNG" alt="Simulink complete model"> 

   The estimation data of the Extended Kalman Filter algorithm for vehicle motion are logged and plotted considering the following scenarios of vehicle motion:

   * Straight motion
   * Curvilinear motion
   * General motion.

   RMSE values of different scenarios using CTRV model: 

   |State|1st scenario RMSE|2nd scenario RMSE| 3rd scenario RMSE|
   |-----------|--------|--------|--------|
   |x|0.6549|0.3949|0.4253|
   |y|0.4893|0.3520|0.3342|
   |theta|0.0149|0.0143|0.0147|
   |v|1.8193|0.8456|0.8304|
   |theta dot|0.0123|0.0168|0.0190|

   Estimation quality may need to be refined (reduce RMSE Error). One possible option is switching to the CTRA motion model.

### Basic Build Instructions
- Install _MATLAB_ (recommended version R2018a),  _DSP System Toolbox_ and _Communications Toolbox_.
- Clone this repo.
- Switch MATLAB current folder (default path) to the project main path.
- Run _Bicycle_CTRV.slx_ or _Bicycle_CTRA.slx_ using _Simulink_ environment.
- Start the simulation using the green (play) button. After the simulation is done, _StopFcn_ callback calls the _CTRV_plotter.m/CTRA_plotter.m_ script to plot and save figures of logged data as PNGs.
- Switch between different scenarios/patterns using the _Signal Builder_ (most left block).
