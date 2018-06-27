% This is just a function that calculates the mean number of function
% evaluations (average) required for an algorithm to converge to the best
% solution as descrived in the DE book. This is the sum of all the function
% evaluation values divided by the number of trials. If visualy compared to
% the maximum number of function evaluations permited it can give a sense
% of the speed of the algorithms. This is obiusly done considering only the
% succesful ones (descard the ones that have not converged)

function [] = MeanFunctEval (label, results, colFE, maxFV)% Inputs: label --> mark or name to be given to the results so it can be uploaded; results --> results from AMIGO2 simulation to be analysed; colFE --> column in the array were the number of function evaluations is located; maxFE --> maximum value of function evaluations that where conceded to all the simulations; 

FE = []; % Empty vector that will include all the values of the function evaluations (all the trials)

for a = 1:length(results) % For loop that goes through the length of the Matrix introduced
    s = results(a,:); % This element extracts each vector (row) of the matrix
    e = s(colFE); % This extracts the selected column for the vector the loop is on
    FE = [FE e]; % Addition of each value of function evaluations at each loop while it goes from vector to vector
end

MFE = (sum(FE)/length(FE)); % Calculation of the mean number of function evaluations as defined
maxFEr = max(FE); % Maximum number of function evaluations used by the algorithm
minFEr = min(FE); % Minimum number of function evaluations used by the algorithm

mkdir(strcat('MeanFunctionEvaluations',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Makes a new  folder that will contain the results
cd(strcat('MeanFunctionEvaluations',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Changes directory to the new folder

save([label,'-MeanFuncEval'],'MFE'); % This saves the value of the average number of fuunction evaluations in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Average number of function evaluations is ', num2str(MFE)]); % This just prints the result on the screan to check it
disp(['The maximum number of function evaluations required by the algorithm was ',num2str(maxFEr)]);
disp(['The minimum number of function evaluations required by the algorithm was ',num2str(minFEr)]);


% prt = (MFE/maxFV)*100; % Percentatge of function evaluations nedded respect the maximum setting selected

y = histogram(FE);
savefig([label, '-histogramMFE.fig']);
label = categorical({'MFE','Maximum Boundary Assigned'}); % Labels for the bar charts
bar(label,[MFE maxFV],'r','FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5); % Bar chart that shows in one bar the average number of function evaluations used and in the other bar the maximum boundary assigned to all the algorithm for a visual representation of how many FE an algorithm needed in average in respect to the total

cd ..\; % Back to the original directory


end