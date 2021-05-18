% Load model file corresponding to a cell of this type
% Has the variables: current, SOC, time, voltage
addpath E:\BMS\SOC\matlab
load E:\BMS\SOC\matlab\PAN_CAPSTONE_DATA.mat; % load data from Panasonic NMC cell, +25 degC
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
J1bestval = 0.15;
J2bestval = 2;
% Get tuning values from user-modified function
for num = 1:200
    [SigmaW, SigmaV, SigmaZ0] = tuneEKF;
    SigmaW = SigmaW(num);
    SigmaV = SigmaV(num);
    SigmaZ0 = SigmaZ0(num);
    SigmaX0 = diag([1e-6 1e-6 SigmaZ0]);
    ekfData = initEKF(voltage(1),T,SigmaX0,SigmaV,SigmaW,model);
    % This simulation tests the EKF when there is an inital SOC-estimation error
    % The true initial SOC is 95%, but we will initialize the SOC estimate in the
    % filter to 90% and see how quickly and well the filter converges toward the
    % correct SOC.
    ekfData.xhat(ekfData.zkInd)=0.90; %
    
    % Now, enter loop for remainder of time, where we update the SPKF
    % once per sample interval
%     fprintf('Please be patient. This code will take a minute or so to execute.\n')
    for k = 1:length(voltage),
        vk = voltage(k); % "measure" voltage
        ik = current(k); % "measure" current
        Tk = T;          % "measure" temperature
        
        % Update SOC (and other model states)
        [sochat(k),socbound(k),ekfData] = iterEKF(vk,ik,Tk,deltat,ekfData);
%         if mod(k,300)==0,
%             fprintf('  Completed %d out of %d iterations...\n',k,length(voltage));
%         end
    end
    
    
    % subplot(1,2,1); plot(time/60,100*sochat,time/60,100*soc); hold on
    % plot([time/60; NaN; time/60],[100*(sochat+socbound); NaN; 100*(sochat-socbound)],'--');
    % title('SOC estimation using EKF'); grid on
    % xlabel('Time (min)'); ylabel('SOC (%)'); legend('Estimate','Truth','Bounds');
    
    
    J1 = sqrt(mean((100*(soc-sochat)).^2));
    % fprintf('RMS SOC estimation error = %g%%\n',J1);
    
    J2 = 100*socbound(end);
    if J1<J1bestval
        J1bestval = J1;
        fprintf('%g,%g,%g,%g\n',SigmaW,SigmaV,SigmaZ0,J1bestval);
    end
    if J2<J2bestval
        J2bestval = J2;
        fprintf('%g,%g,%g,%g\n',SigmaW,SigmaV,SigmaZ0,J2bestval);
    end
end
%%
% fprintf('Final value of SOC estimation error bounds = %g%%\n',J2);

% subplot(1,2,2); plot(time/60,100*(soc-sochat)); hold on
% plot([time/60; NaN; time/60],[100*socbound; NaN; -100*socbound],'--');
% title('SOC estimation errors using EKF');
% xlabel('Time (min)'); ylabel('SOC error (%)'); ylim([-4 4]); 
% legend('Estimation error','Bounds'); 
% grid on

ind = find(abs(soc-sochat)>socbound);
fprintf('Percent of time error outside bounds = %g%%\n',length(ind)/length(soc)*100);

% Compute the prospective grade
tableRow = min(11,ceil(max(0,J1-0.1)/0.01 + 1));
tableCol = min(11,ceil(max(0,J2-0.21)/0.02 + 1));
table = hankel([10:-1:0]);
grade = table(tableRow,tableCol);
if ~isempty(ind),
  fprintf('Your SOC estimation error was sometimes outside of bounds, so your overall grade is 0/10.');
else
  fprintf('Your grade is calculated from row %d and column %d of the grading table that is\n',tableRow,tableCol);
  fprintf('listed in the project description page. This will result in a grade of %d/10.\n',grade);
end