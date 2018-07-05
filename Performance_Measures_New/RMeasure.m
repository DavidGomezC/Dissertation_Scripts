function [] = RMeasure (measure, label)% Inputs: label --> mark or name to be given to the results so it can be uploaded; results --> results from AMIGO2 simulation to be analysed; best --> The value (Qm, p, C, mean cost function) of the tunning with the best results

%%%%%%%%%%%%%%
%  DISTANCE  %
%%%%%%%%%%%%%%

cd (measure);
% load('test-Q_Measure');
SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),label)) 
        list2=[list2; SN(x,:)]; 
    end
end

vper = [];
best = [];

for x=1:length(list2(:,1))
    load(list2(x,:));
    if strcmp(measure, 'Q_Measure')
        vper = [vper, Qm];
        best = max(vper);
    elseif strcmp(measure, 'AverageOfConvergence')
        vper = [vper, p];
        best = max(vper);    
    elseif strcmp(measure, 'MeanFunctionEvaluations')
        vper = [var, MFE];
        best = min(vper);
    end
end

cd ..\;

mn = sum(vper)/length(vper); % Mean number of the results from the performance measure used to calculate the robustness
dis=abs(best-mn); % Distance from the best value to the mean of values
disp(['Best value is: ',num2str(best)]); % Printing of the results on the screan
disp(['The mean value is: ',num2str(mn)]);
disp(['The difference between the two values (distance) is: ',num2str(dis)]);
fprintf('\n*\nNEXT\n*\n');
 
%%%%%%%%%%%%%%%%%%%%%%%%
%  STANDARD DEVIATION  %
%%%%%%%%%%%%%%%%%%%%%%%%

par=[]; % Empty vector
for c=1:length(vper) % Goes through the length of the vector results
    e=(vper(c)-mn)^2; % Substracts mean and squares the result for each element in the results
    par = [par e]; % Adds each result into the empty vector
end
sd = sqrt(sum(par)/length(par)); % Calculation of the standard deviation

cd(strcat('R-Measure')); % Changes directory to the new folder

save([label,'-Standard Deviation'],'sd'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['The standard deviation is ', num2str(sd)]); % This just prints the result on the screan to check it
fprintf('\n*\nNEXT\n*\n');

%%%%%%%%%%%%%%
%  VARIANCE  %
%%%%%%%%%%%%%%

vr = (sd)^2; % Calculation of the Variance

save([label,'-Variance'],'vr'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['The variance is ', num2str(vr)]); % This just prints the result on the screan to check it
fprintf('\n*\nNEXT\n*\n');

%%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFIED DEVIATION  %
%%%%%%%%%%%%%%%%%%%%%%%%

parm=[]; % Empty vector
for c=1:length(vper) % Goes through the length of the vector results
    e=(vper(c)-best)^2; % Substracts mean and squares the result for each element in the results
    parm = [parm e]; % Adds each result into the empty vector
end
msd = sqrt(sum(parm)/length(parm)); % Calculation of the standard deviation

save([label,'-ModifiedStandardDeviation'],'msd'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['The modified standard deviation is ', num2str(msd)]); % This just prints the result on the screan to check it
fprintf('\n*\nNEXT\n*\n');

%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFIED VARIANCE  %
%%%%%%%%%%%%%%%%%%%%%%%

mvr = (msd)^2; % Calculation of the Variance

save([label,'-ModifiedVariance'],'mvr'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['The modified variiancen is ', num2str(mvr)]); % This just prints the result on the screan to check it

cd ..\; % Back to the original directory
end