% Load model file corresponding to a cell of this type
% Has the variables: current, SOC, time, voltage
addpath E:\BMS\SOC\matlab
load E:\BMS\SOC\matlab\PANdata_P25.mat; % load data from Panasonic NMC cell, +25 degC
T = 25; % Test temperature

time    = DYNData.script1.time(:);   deltat = time(2)-time(1);
time    = time-time(1); % start time at 0
current = DYNData.script1.current(:); % discharge > 0; charge < 0.
voltage = DYNData.script1.voltage(:);
soc     = DYNData.script1.soc(:);

% Load cell-test data to be used for this batch experiment
% Contains variable "DYNData" of which the field "script1" is of 
% interest. This has sub-fields time, current, voltage, soc.
load E:\BMS\SOC\matlab\PANmodel.mat; % load ESC model of Panasonic NMC cell

% Reserve storage for computed results, for plotting
sochat = zeros(size(soc));
socbound = zeros(size(soc));

% Covariance values
SigmaX0 = diag([1e2 1e-2 1e-3]); % uncertainty of initial state
SigmaV = 3e-1; % Uncertainty of voltage sensor, output equation
SigmaW = 4e0; % Uncertainty of current sensor, state equation

% Create spkfData structure and initialize variables using first
% voltage measurement and first temperature measurement
spkfData = initSPKF(voltage(1),T,SigmaX0,SigmaV,SigmaW,model);

% Now, enter loop for remainder of time, where we update the SPKF
% once per sample interval
fprintf('Please be patient. This code will take several minutes to execute.\n')
for k = 1:length(voltage),
  vk = voltage(k); % "measure" voltage
  ik = current(k); % "measure" current
  Tk = T;          % "measure" temperature
  
  % Update SOC (and other model states)
  [sochat(k),socbound(k),spkfData] = iterSPKF(vk,ik,Tk,deltat,spkfData);
  % update waitbar periodically, but not too often (slow procedure)
  if mod(k,1000)==0,
    fprintf('  Completed %d out of %d iterations...\n',k,length(voltage));
  end  
end

%%
subplot(1,2,1); plot(time/60,100*sochat,time/60,100*soc); hold on
plot([time/60; NaN; time/60],[100*(sochat+socbound); NaN; 100*(sochat-socbound)]);
title('SOC estimation using SPKF'); grid on
xlabel('Time (min)'); ylabel('SOC (%)'); legend('Estimate','Truth','Bounds');

%%
fprintf('RMS SOC estimation error = %g%%\n',sqrt(mean((100*(soc-sochat)).^2)));

%%
subplot(1,2,2); plot(time/60,100*(soc-sochat)); hold on
plot([time/60; NaN; time/60],[100*socbound; NaN; -100*socbound]);
title('SOC estimation errors using SPKF');
xlabel('Time (min)'); ylabel('SOC error (%)'); ylim([-4 4]); 
legend('Estimation error','Bounds'); 
grid on

ind = find(abs(soc-sochat)>socbound);
fprintf('Percent of time error outside bounds = %g%%\n',length(ind)/length(soc)*100);