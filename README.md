# BatteryManagementSystem_Algorithm_usingMatlab
**Overview of the code**
 The codes are used in supplement with the books **Battery Management Systems I and II, by Gregory Plett** and the course **Algorithms for Battery Management System, by Coursera**
## State_Of_Charge
 
  **Abbreviations\Symbols used in the code and their physical\actual meaning and units**
  - z    :  State of Charge
  - R    :  Internal Resistance(Ohms)
  - Q    :  Cell capacity(Wh)
  - v    :  Cell  terminal voltage(V)
  - rc   :  Resistor Capacitor Parameters
  - s    :  Instantaneous Hysteresis
  - h    :  Dynamic Hysteresis
  - t    :  Time
 
 
 
 This branch(State_of_charge) of the main folder(BatteryManagementSystem_Algorithm_usingMatlab) aims to determine the State of Charge of an individual Lithium-Ion Cell and also the State of Charge of an Electric Vehicle Battery. The methods used consist of Kalman Filter methods and its varieties. The Kalman filter based approach is the best method for linear systems. While for the nonlinear systems Unscented Kalman Filter based methods are used. To use these algorithms, a good knowledge of Kalman filter theory is required.

**A brief introduction to Kalman filter methods**
 A dynamic physical system, in addition to the input and output signals also is superimposed of noises. The noises can generally be classified as process noise and sensor noise.The process noise are the noise because of the system modelling error,while the sensor noise is the noise picked up by the output measurements.The Kalman filter is a robust method used to reduce the  noise to a minimal amount, which resulting in increasing the accuracy of the system. The Kalman filter generally proceeds in two steps:
- The first step is the **Prediction Step**. In this step the prediction of the states of the current states of the system based on the previous states takes place.
- The second step is the **Correct/Updation Step**. In this the step the future value is estimated with the help of a particular ***Gain Matrix***.
The Gain Matrix is a parameter that disinguishes the Kalman Filter with all the other methods.

**Extended Kalman filter**
In this method,instead of passing all the input points,along with the noise, through the system and obtaining an output, a set of weighted points that accurately represent the behaviour of the system is passed through the system and an output in terms of the weighted points are obtained. The final state of the system is then extrapolated from these points. This method is very useful for nonlinear physical systems.

**Bar Delta Filtering**
One importamt assumption in Kalman filter approach is that all the noises are assumed to be zero mean. But physically one cannot predict the behaviour of noise and generally the noise is not zero mean and there is some current-sensor bias which can cause permanent SOC error.
Also,Accumulated ampere-hours of bias tend to move SOC estimate faster than measurement updates can correct. 
So, to counter these obstacles **Bar Delta Filtering** method is used.
This method estimates the bias and uses this as the input signal to the Kalman Filters.
   
 **Explanation of the codes used**
 
 ***Extended explanation of the code is present in the comments of the code file***
 1) Code names: **KFinaction** **EKFinAction.m** 
    - **Explanation** This code is used to simulate the Extended Kalman filter on a simple system.
 
 2) Code name: **initEKF.m** **initSPKF.m** **initSPKFbd**
    - **Explanation** This code initializes interprets the data and initializes the variable.
 
 3) Code name : **iterEKF.m** **iterSPKF.m** **iterSPKFbd**
    - **Explanation** This is the section where the actual Filter process takes place and the output is stored.
 
 4) Code name : **wrapperEKF.m** **wrapperSPKF.m** **wrapperSPKFbd**
    - **Explanation** This section simulates the entire process and plots and visualizes the data.
 
 5) Code name : **simCell.m**
    - **Explanation** This  simulates the behaviour of a single cell.

 6) Code name : **bardeltadata.m**
    - **Explanation** This provides  with the data for the Bar delta filtering method.
   
## ElectroChemical Cell Model

This branch of algorithm aims to exexute the following things.
- Interpret the data obtained from the laboratory tests on a lithium-ion cell, and use this data to obtain the cell parameters,for example: State of Charge, Terminal voltage, Open Circuit Voltage for a cell.
The cell is modelled as a Randless circuit with Warburg impedance.
![Image of Randless circuit ](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/randlesscircuit.JPG)
The equations for this model are formulated, first in continuous time model and then converted to discrete time model.
The currents from the resistor current pairs are 
![Image of RC current ](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/resistor_current.JPG)
**Hysteresis** is an important parameter that needs to be taken into account while modelling a lithium-ion cell. In an normal cell, if it is allowed to rest long enough, diffusion voltagesdecay to zero, so model voltage decays to OCV, but due to hysteresis, the behaviour of lithium-ion cell is different.
The hysteresis is differentiated as two forms and consequently, the equations are formulated
- a) Dynamic hysteresis 
- ![Image of dynamic hysteresis](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/hysteresis.JPG)
- b) Instantaneous hysteresis 
     - ![Image of instantaneous hysteresis](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/instantaneoushysteresis.JPG)
       The **instantenous hysteresis** only changes when there is a change in the direction of the current.
        The overall hysteresis is the sum of the two hysteresis effects.
 
 If all the effects are cumulated, the cell model is said to be an **Enhanced Self Correcting Model**.
 The state equation of an ESC is 
 ![Image of an ESC state equation](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/modelequation.JPG)
 
 The ESC output equation is 
 ![Image of ESC output equation](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/opequation.JPG)
 
 Finally, the state space form of an ESC model is 
 ![Image of State space form](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/statespaceformofcell.JPG).
 
 Code **runProcessOCV.m** runs the processOCV for a cell against the recorded laboratory data.
 
 While **processOCV** simulates the cell for static data/cell without a load, the **processDynamicOCV** simulates the cell for dynamic data/cell under load.
 
 The overall procedure to deteremine parameter values for a cell under load is as follows:
 1. First, ***Coulombic Efficiency*** and ***Total Cell Capacity*** are calculated from data, directly, as we did for the OCV test results
 2. Compute ***state of charge*** and ***open circuit voltage*** for every data sample; subtract OCV from terminal voltage.
 3. Use ***subspace system identiﬁcation*** technique to ﬁnd R–C time constants.
 4. Compute ***instantaneous hysteresis*** and ***R-C currents*** for every data sample
 5. A ***hystereis rate constant*** is calculated that determines the amount of hysteresis present in a cell.
 6. "Unexplained” part of voltage is now linear in parameters—solve for these parameter values using least squares.
    - The unexplained part is
    - ![Image of unexplained part](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/unexplainedpart.JPG)
    - ![Image of lsq method matlab](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/unexplainedpart2.JPG)
    - The X matrix is calculated in least square method.
 7. Compute rms voltage-prediction error of present model
 8. Adapt ***hysteresis rate constant*** to minimize this error, iterating steps 5–8 until convergence is reached

 Code **getParamESC.m** retrieves a parameter for a given cell at a given temperature.
 
 Code **OCVfromSOCtemp.m** interpolates and returns the open circuit voltage for a given state of charge and temperature.
 
 Simulation of a Constant current/Constant Voltage charging cycle can be performed with **ChargingCycleSimulation.m**
 
 Simulation of a Constant power/Constant Voltage charging cycle can be performed with **ChargingCycleSimulation(CP_CV).m**
 
 Simulation of a Parallel Cell module can be performed with **PCMsim.m**
 
 Simulation of a Series Cell module can be performed with **SCMsim.m**
 
 To take the algorithms a step further, the battery is simulated for an electric vehicle under various drive cycles and taking under consideration the characteristics of motor,inverter,drivetrain,etc.
The drive cycles used are **Urban Dynamometer Driving Schedule**, **Highway Fuel Efficiency Test**, **New York City Cycle**.

The approach to simulating an electric vehicle behaviour is as follows:

![Image of EVsimulation](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/blockdiagram.JPG)

  - The desired accelerations and desired speed are limited by the actual achievable accelerations and speed of the vehicle and powertrain module.
  - The desired road forces are limited by the forces that the tire-road dynamic can actually achieve.
  - Desired torque and power values restricted by speciﬁcations of motor chosen for vehicle, and thus achievable torques and power values may be lower. The desired  torque must be compatible with the speed-torque for the particular motor chosen.
  - Consequently the the battery power is computed based on these limitations and the limitation of the motor power.
  - The Battery SOC is then updated.
  - Finally the Battery range is extrapolated from the rate of depletion of the availabe battery energy.
 
 **EQUATIONS USED**:
 
   - The vehicle desired acceleration is
     ![Image of desired acceleration](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/desireedacc.JPG)
     
   - The desired acceleration force is
     ![Image of desired acceleration force](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/desaccforce.JPG)
       - The equivalent mass required is 
         ![Image of equivalent mass](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/eqmass.JPG)
         
   - The drag forces acting on the vehicle are classified as aerodynamic force and rolling force 
     ![Image of drag forces](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/dragforces.JPG)
     
   - The brake drag and grade forces acting on the vehicle are
     ![Image of brake drag and rolling force](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/dragf2.JPG)
     
   - The demanded motor torque can be derived as 
     ![Image of demanded motor torque](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/desireedacc.JPG)
     
   - The actual acceleration force and consequently the actual acceleration is as follows
     ![Image of actual acceleration and force](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/accforce.JPG)
     
   - The limit motor speed is 
     ![Image of limit motor speed](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/motorspeed.JPG)
 
   - The actual vehicle vehicle speed is then computed from this limit speed
     ![Image of vehicle speed](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/actualspeed.JPG)
     
   - The battery power demand is then computed
     ![Image of motor power](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/mototrpower.JPG)
     ![Image of battery power](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/battpower.JPG)
     
   - The current drawn out of the battery is 
     ![Image of battery current](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/battcurrent.JPG)
     
   - The State of Charge of the battery for the drive cycle is 
     ![Image of battery SOC](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/battSOC.JPG)
     
   - Finally the driving range of the battery according to the different drive cycles is calculated
     ![Image of battery range](https://github.com/Ronakj2904/BatteryManagementSystem_Algorithm_usingMatlab/blob/master/images/range.JPG)
     
 Code **setupSimVehicle.m** sets up the different blocks for the simulation of the electric vehicle behaviour. A different block is used each for setting up a single cell,PCM/SCM module, Battery pack, Motor, Wheel, Drivetrain and Vehicle.
 
 Code **simVehicle.m** simulates the Vehicle for the given data and stores the result.
 
 Code **setupSimVehicleplots.m** plots the obtained results.
 
 
