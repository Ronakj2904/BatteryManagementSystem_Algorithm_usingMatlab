%%
% Plot the desired and actual vehicle speed profiles
names = {'NYCC','UDDS','US06','HWY'};
style = {'b','r','m','g'};
locations = {'northwest','northeast','south','south'};
for ii = 1:length(results)
  plotData = results{ii}.actualSpeedKPH;
  subplot(2,2,ii); plot(0:length(plotData)-1,plotData,style{ii});
  plotData = results{ii}.desSpeedKPH; hold on; grid on;
  plot(0:length(plotData)-1,plotData,'k');
  title(sprintf('%s profile: vehicle speed, 0.3%% grade',names{ii}));
  xlabel('Time (s)'); ylabel('Speed (km h^{-1})');
  legend('Actual','Desired','location',locations{ii});
end

%%
% Plot the Battery Pack Current Profiles
for ii = 1:length(results),
  plotData = results{ii}.current;
  subplot(2,2,ii); plot(0:length(plotData)-1,plotData,style{ii}); 
  title(sprintf('%s profile: battery current',names{ii})); 
  xlabel('Time (s)'); ylabel('Current (A)'); grid on;
end
%% Plot the Battery Pack Demanded Power Profile
for ii = 1:length(results),
  plotData = results{ii}.batteryDemand;
  subplot(2,2,ii); plot(0:length(plotData)-1,plotData,style{ii}); 
  title(sprintf('%s profile: battery power',names{ii})); 
  xlabel('Time (s)'); ylabel('Power (kW)'); grid on;
end
%% Plot Histograms of Battery Pack Demanded Power
edges = -70:20:140;
for ii = 1:length(results),
  [N,Bin] = histc(results{ii}.batteryDemand,edges);
  subplot(2,2,ii); bar(edges,N,'histc'); xlim([min(edges) max(edges)]);
  title(sprintf('%s profile: battery power',names{ii}));
  ylabel('Frequency'); xlabel('Power (kW)'); grid on;
end
%% Plot Histograms of Motor Demanded Power
for ii = 1:length(results)
  [N,Bin] = histc(results{ii}.limitPower,edges);
  subplot(2,2,ii); bar(edges,N,'histc'); xlim([min(edges) max(edges)]);
  title(sprintf('%s profile: limited motor power',names{ii}));
  ylabel('Frequency'); xlabel('Power (kW)'); grid on;
end
%% Plot Scatter Plot of Motor Torque Vs. RPM
motor = results{1}.vehicle.drivetrain.motor;
regen = results{1}.vehicle.drivetrain.regenTorque;
RPM = motor.RPMrated:50:motor.RPMmax;
TRQ = motor.Lmax*motor.RPMrated./RPM;
for ii = 1:length(results),
  subplot(2,2,ii);
  scatter(results{ii}.motorSpeed,results{ii}.limitTorque,style{ii});
  title(sprintf('%s profile: motor torque vs. speed',names{ii}));
  xlim([0 12000]); ylim([-300 300]); hold on; grid on;
  xlabel('Speed (RPM)'); ylabel('Torque (N m)');
  plot([0 motor.RPMrated],motor.Lmax*[1 1],'k');
  plot([0 motor.RPMrated],-regen*motor.Lmax*[1 1],'k');
  plot(RPM,TRQ,'k'); plot(RPM,-regen*TRQ,'k');
end