


% function [SigmaW, SigmaV, SigmaZ0] = tuneEKF
%
% SigmaW - covariance value for current-sensor process noise
% SigmaV - covariance value for voltage-sensor measurement noise
% SigmaZ0 - covariance value for error in initial SOC estimate

function [SigmaW, SigmaV, SigmaZ0] = tuneEKF

  % BEGIN MODIFYING CODE AFTER THIS
  SigmaW = linspace(0.0005,0.0008,200);
  SigmaV = linspace(0.0006,0.0009,200);
  SigmaZ0 = linspace(0.01,0.015,200);
 
  end
