**Overview of the code**
 The codes are used in supplement with the books **Battery Management Systems I and II, by Gregory Plett** and the course **Algorithms for Battery Management System, by Coursera**
# State_Of_Charge
 
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
   
# ElectroChemical Cell Model

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


