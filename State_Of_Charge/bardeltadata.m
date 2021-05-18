% Use a desktop-validation approach to create a synthetic dataset
% This dataset assumes a 4-cell battery pack comprising Panasonic 25 Ah
% cells. The cells have different series resistances, different capacities,
% and different initial states of charge. They are all exercised with the
% same profile of current versus time, which comprises two UDDS profiles
% separated by rest intervals.
addpath E:\BMS\SOC\matlab                  % set your own path to the algorithm/code
load E:\BMS\SOC\matlab\PANmodel;           % ESC model of Panasonic cell
udds = load('E:\BMS\SOC\matlab\udds.txt'); % profile of current versus time for UDDS
T = 25;

z0 = 0.9:-0.1:0.6;       % set initial SOC for each cell in pack
R0 = (1.3:-0.1:1)*1e-3;  % set R0 for each cell in pack
Q0 = 25:28;              % set Q for each cell in pack
% create current profile: rest/udds/rest/udds/rest
ik = [zeros(300,1); udds(:,2); zeros(300,1); udds(:,2); zeros(241,1)];

vk = zeros(length(ik),length(z0)); % reserve storage for cell voltages
zk = zeros(length(ik),length(z0)); % reserve storage for cell SOCs
    
for k = 1:length(z0),
  model.R0Param = R0(k)*ones(size(model.R0Param)); % overwrite R0
  model.QParam  = Q0(k)*ones(size(model.QParam)); % overwrite Q
  [vcell,rck,hk,zcell,sik,OCV] = simCell(ik,T,1,model,z0(k),0,0);
  vk(:,k) = vcell;
  zk(:,k) = zcell;
end
save('PANPackData.mat','vk','zk','ik','T','Q0','R0');

subplot(1,2,1); t = (0:length(ik)-1)/3600;
plot(t,vk); xlabel('Time (hr)'); ylabel('Voltage (V)');
title('Voltage versus time for four cells');

subplot(1,2,2);
plot(t,100*zk); xlabel('Time (hr)'); ylabel('State of charge (%)');
title('SOC versus time for four cells');