% It measures the robustness of an algorythm and it is good to evaluate
% different tunings of an algorithm. NONE of the two ways to calculate it 
% can be used since it is based in comparison of the best tuning parameters
% between different benchmark problem tests and we are only dealing with
% one problem here. A different way to measure robustnes of an algorithm is
% introduced here then:

% Hypothesis --> I am going to use instead of the tunning paramters the
% results from each tuning, and for the measure either Q-measure,
% probability of convergence or convergence measure can be used, this
% depend on what you want to compare. The mean value of the cost function
% can be used as well. The tunning with the best result has to be
% identified a priory. With this the function will calculate the mean value
% and compare it to the best value (to see distance from them), standard
% deviation (how spread out are the values respect the mean), variance (how
% far are the values from the mean) but also calculate standard deviation
% and variance in respect from the best value (not the mean) where if the
% algorithm is robust (different tunnings will lead to the same result)
% this last two will be close to 0 but if it is not robust (different 
% tunnings lead to different results) then this two will have a larger
% value (we want the minimum value possible)

function [] = RMeasure (label, results, best)% Inputs: label --> mark or name to be given to the results so it can be uploaded; results --> results from AMIGO2 simulation to be analysed; best --> The value (Qm, p, C, mean cost function) of the tunning with the best results

%%%%%%%%%%%%%%
%  DISTANCE  %
%%%%%%%%%%%%%%
mn = sum(results)/length(results); % Mean number of the results from the performance measure used to calculate the robustness
dis=abs(best-mn); % Distance from the best value to the mean of values
disp(['Best value is: ',num2str(best)]); % Printing of the results on the screan
disp(['The mean value is: ',num2str(mn)]);
disp(['The difference between the two values (distance) is: ',num2str(dis)]);


%%%%%%%%%%%%%%%%%%%%%%%%
%  STANDARD DEVIATION  %
%%%%%%%%%%%%%%%%%%%%%%%%
% I had some isues with the matlab function std and var so I will do it
% manualy

par=[]; % Empty vector
for c=1:length(results) % Goes through the length of the vector results
    e=(results(c)-mn)^2; % Substracts mean and squares the result for each element in the results
    par = [par e]; % Adds each result into the empty vector
end
sd = sqrt(sum(par)/length(par)); % Calculation of the standard deviation

mkdir(strcat('R-Measure',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Makes a new  folder that will contain the results
cd(strcat('R-Measure',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Changes directory to the new folder

save([label,'-Standard Deviation'],'sd'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['The standard deviation is ', num2str(sd)]); % This just prints the result on the screan to check it


%%%%%%%%%%%%%%
%  VARIANCE  %
%%%%%%%%%%%%%%

vr = (sd)^2; % Calculation of the Variance

save([label,'-Variance'],'vr'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['The variiancen is ', num2str(vr)]); % This just prints the result on the screan to check it


%%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFIED DEVIATION  %
%%%%%%%%%%%%%%%%%%%%%%%%

parm=[]; % Empty vector
for c=1:length(results) % Goes through the length of the vector results
    e=(results(c)-best)^2; % Substracts mean and squares the result for each element in the results
    parm = [parm e]; % Adds each result into the empty vector
end
msd = sqrt(sum(parm)/length(parm)); % Calculation of the standard deviation

save([label,'-ModifiedStandardDeviation'],'msd'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['The modified standard deviation is ', num2str(msd)]); % This just prints the result on the screan to check it

%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFIED VARIANCE  %
%%%%%%%%%%%%%%%%%%%%%%%

mvr = (msd)^2; % Calculation of the Variance

save([label,'-ModifiedVariance'],'mvr'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['The modified variiancen is ', num2str(mvr)]); % This just prints the result on the screan to check it



cd ..\; % Back to the original directory
end