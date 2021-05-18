% results = simVehicle(vehicle,cycle,grade)
%   - simulate a vehicle defined by "vehicle", perhaps created using
%     setupSimVehicle.m
%   - cycle is Nx2, where first column is time in seconds and second column
%     is desired speed in miles per hour
%   - grade is road grade in percent - either a constant grade for all
%     time, or a different grade value for every point in time
function results = simVehicle(vehicle,cycle,grade)
  rho = 1.225; % air density, kg/m3

  results.vehicle = vehicle;
  results.cycle = cycle; % time in s, desired speed in miles/hour
  results.time = cycle(:,1); % s
  if isscalar(grade),
    results.grade = repmat(atan(grade/100),size(results.time)); % rad
  else
    results.grade = atan(grade/100); % rad
  end 
  results.desSpeedKPH = cycle(:,2) * 1.609344; % convert to km/h
  results.desSpeed = min(vehicle.maxSpeed,results.desSpeedKPH*1000/3600); % m/s

  % pre-allocate storage for results
  results.desAccel = zeros(size(results.desSpeed)); % m/s2
  results.desAccelForce = zeros(size(results.desSpeed)); % N
  results.aeroForce = zeros(size(results.desSpeed)); % N
  results.rollGradeForce = zeros(size(results.desSpeed)); % N
  results.demandTorque = zeros(size(results.desSpeed)); % N-m
  results.maxTorque = zeros(size(results.desSpeed)); % N-m
  results.limitRegen = zeros(size(results.desSpeed)); % N-m
  results.limitTorque = zeros(size(results.desSpeed)); % N-m
  results.motorTorque = zeros(size(results.desSpeed)); % N-m
  results.demandPower = zeros(size(results.desSpeed)); % kW
  results.limitPower = zeros(size(results.desSpeed)); % kW
  results.batteryDemand = zeros(size(results.desSpeed)); % kW
  results.current = zeros(size(results.desSpeed)); % A
  results.batterySOC = zeros(size(results.desSpeed)); % 0..100
  results.actualAccelForce = zeros(size(results.desSpeed)); % N
  results.actualAccel = zeros(size(results.desSpeed)); % m/s2
  results.motorSpeed = zeros(size(results.desSpeed)); % RPM
  results.actualSpeed = zeros(size(results.desSpeed)); % m/s
  results.actualSpeedKPH = zeros(size(results.desSpeed)); % km/h
  results.distance = zeros(size(results.desSpeed)); % km
  
  prevSpeed = 0; prevTime = results.time(1) - 1; prevMotorSpeed = 0; 
  prevSOC = vehicle.drivetrain.pack.socFull; prevDistance = 0;
  for k = 1:length(results.desSpeed),
    results.desAccel(k) = (results.desSpeed(k) - prevSpeed)/ ...
        (results.time(k) - prevTime);
    results.desAccelForce(k) = vehicle.equivMass * results.desAccel(k);
    results.aeroForce(k) = 0.5 * rho * vehicle.Cd * vehicle.A * prevSpeed^2;
    results.rollGradeForce(k) = vehicle.maxWeight * 9.81 * sin(results.grade(k));
    if abs(prevSpeed) > 0,
      results.rollGradeForce(k) = results.rollGradeForce(k) + ...
        vehicle.drivetrain.wheel.rollCoef * vehicle.maxWeight * 9.81;
    end
    results.demandTorque(k) = (results.desAccelForce(k) + results.aeroForce(k) + ...
        results.rollGradeForce(k) + vehicle.roadForce) * ...
        vehicle.drivetrain.wheel.radius / vehicle.drivetrain.gearRatio;
    if prevMotorSpeed < vehicle.drivetrain.motor.RPMrated,
      results.maxTorque(k) = vehicle.drivetrain.motor.Lmax;
    else
      results.maxTorque(k) = vehicle.drivetrain.motor.Lmax * ...
          vehicle.drivetrain.motor.RPMrated / prevMotorSpeed;
    end
    results.limitRegen(k) = min(results.maxTorque(k),...
        vehicle.drivetrain.regenTorque  * vehicle.drivetrain.motor.Lmax);
    results.limitTorque(k) = min(results.demandTorque(k),results.maxTorque(k));
    if results.limitTorque(k) > 0,
      results.motorTorque(k) = results.limitTorque(k);
    else
      results.motorTorque(k) = max(-results.limitRegen(k),results.limitTorque(k));
    end
    
    results.actualAccelForce(k) = results.limitTorque(k) * vehicle.drivetrain.gearRatio / ...
        vehicle.drivetrain.wheel.radius - results.aeroForce(k) - results.rollGradeForce(k) - ...
        vehicle.roadForce;
    results.actualAccel(k) = results.actualAccelForce(k) / vehicle.equivMass;
    results.motorSpeed(k) = min(vehicle.drivetrain.motor.RPMmax,...
        vehicle.drivetrain.gearRatio * (prevSpeed + results.actualAccel(k) * ...
        (results.time(k) - prevTime)) * 60 / (2*pi*vehicle.drivetrain.wheel.radius));
    results.actualSpeed(k) = results.motorSpeed(k) * 2*pi*vehicle.drivetrain.wheel.radius / ...
        (60 * vehicle.drivetrain.gearRatio);
    results.actualSpeedKPH(k) = results.actualSpeed(k) * 3600/1000;
    deltadistance = (results.actualSpeed(k) + prevSpeed)/2 * (results.time(k) - prevTime)/1000;
    results.distance(k) = prevDistance + deltadistance;

    if results.limitTorque(k) > 0,
      results.demandPower(k) = results.limitTorque(k);
    else
      results.demandPower(k) = max(results.limitTorque(k),-results.limitRegen(k));
    end
    results.demandPower(k) = results.demandPower(k) * 2*pi * (prevMotorSpeed + results.motorSpeed(k))/2 / 60000;
    results.limitPower(k) = max(-vehicle.drivetrain.motor.maxPower,min(...
        vehicle.drivetrain.motor.maxPower,results.demandPower(k)));
    results.batteryDemand(k) = vehicle.overheadPwr/1000;
    if results.limitPower(k) > 0,
      results.batteryDemand(k) = results.batteryDemand(k) + ...
            results.limitPower(k)/vehicle.drivetrain.efficiency;
    else
      results.batteryDemand(k) = results.batteryDemand(k) + ...
            results.limitPower(k)*vehicle.drivetrain.efficiency;
    end
    results.current(k) = results.batteryDemand(k)*1000/vehicle.drivetrain.pack.vnom;
    results.batterySOC(k) = prevSOC - results.current(k) * (results.time(k) - prevTime) / ...
        (36*vehicle.drivetrain.pack.module.capacity);
      
    
    prevTime = results.time(k);
    prevSpeed = results.actualSpeed(k);
    prevMotorSpeed = results.motorSpeed(k);
    prevSOC = results.batterySOC(k);
    prevDistance = results.distance(k);
  end
end
