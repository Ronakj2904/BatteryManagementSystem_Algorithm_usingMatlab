clear
% Can set an artificial current-sensor bias to see how well the SPKF
% estimates this bias
ibias = 0.5; % can set this (e.g., to 0 or 0.5 A)

% Load cell-test data . Contains variable "DYNData" of which the field
% "script1" is of interest. It has sub-fields time, current, voltage, soc.
addpath E:\BMS\SOC\matlab
load E:\BMS\SOC\matlab\PANmodel.mat; % loads cell model
load E:\BMS\SOC\matlab\PANPackData.mat;
time = 0:length(ik)-1; time = time(:); deltat = 1;
current = ik + ibias;
voltage = vk;
soc = zk;

% For this somewhat simplified example, we will not estimate 
% capacity inverse. Instead, we assume perfect knowledge of capacity
% of every cell, and hence capacity inverse of every cell
Qinv = 1./Q0; Qinvbar = mean(Qinv); dQinv = Qinv - Qinvbar;
% Similarly, we will not estimate the delta-R0 values, but assume knowledge
dR0 = getParamESC('R0Param',T,model) - R0;

% Reserve storage for computed results, for plotting
sochat = 0*time;    % reserve storage for bar-soc values
socbound = sochat;  % and bounds on those values
bias = sochat;      % ... also for current-sensor bias estimate
biasBound = sochat; % and bounds on that estimate

dsochat = zeros(size(voltage));   % reserve storage for delta-soc values
dsocbound = zeros(size(voltage)); % and for bounds on those values
dR = zeros(size(voltage));        % reserve storage for delta-R0 values
dRbound = zeros(size(voltage));   % and for bounds on those values

% Covariance values
% State ordering: ir,h,z,bias
SigmaX0 = diag([1e2 1e-4 1e-2 5e-2]); % uncertainty of initial state
SigmaV = 1e-3; % uncertainty of voltage sensor, output equation
SigmaW = diag([1e-1, 1e-4]); % uncertainty of current sensor, bias

spkfData = initSPKFbd(voltage(1,:), T, SigmaX0, SigmaV, SigmaW, model);
spkfData.Qinvbar = Qinvbar; spkfData.dQinv = dQinv; spkfData.dR0 = dR0;

% Now, enter loop for remainder of time, where we update the SPKF
% once per sample interval
fprintf('Starting SPKF\n');
for k = 1:size(voltage,1),
  vk = voltage(k,:); % "measure" voltage
  ik = current(k); % "measure" current
  Tk = T; % "measure" temperature

  % Update SOC (and other model states) of bar filter
  [sochat(k), socbound(k), dsochat(k,:), dsocbound(k,:), bias(k), biasBound(k), spkfData] = ...
      iterSPKFbd(vk, ik, Tk, deltat, spkfData);

  % update waitbar periodically, but not too often (slow procedure)
  if mod(k,250)==0, fprintf('  Completed %d out of %d iterations\n',k,size(voltage,1)); end
end

% Display output
subplot(2,2,1); 
plot(time/60,100*sochat, '-r',  time/60, 100*mean(soc,2), '--b'); hold on
h = plot([time/60; NaN; time/60],...
        [100*(sochat+socbound); NaN; 100*(sochat-socbound)]);
plot(time/60,100*soc, '-k');
ylim([55 95]);
title('Avg. SOC estimate using bar filter' );
xlabel('Time (min)' ); ylabel('SOC (%)');
legend('Average SOC estimate', 'True average SOC', 'Bounds', 'Individual SOCs'); grid on

subplot(2,2,2);
plot(time/60,bias, '-r',  time/60, ibias*ones(size(bias)), '--b'); hold on
h = plot([time/60; NaN; time/60],...
        [bias+biasBound; NaN; bias-biasBound]);
title('Current-sensor bias estimate using bar filter' );
xlabel('Time (min)' ); ylabel('Bias (A)' );
legend('Bias estimate', 'True bias', 'Bounds'); grid on

subplot(2,2,3);
sochat = repmat(sochat,1,4);
socbound = repmat(socbound,1,4);
plot(time/60,100*(sochat+dsochat), '-r'); hold on
h = plot([time/60; NaN; time/60],...
        [100*(sochat+dsochat+socbound+dsocbound); NaN(1,4); 100*(sochat+dsochat-dsocbound-socbound)]);
plot(time/60,100*soc, '-k');
ylim([55 95]);
title('Individual SOC estimate using bar-delta filter' );
xlabel('Time (min)' ); ylabel('SOC (%)'); grid on