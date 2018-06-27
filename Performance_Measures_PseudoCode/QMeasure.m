% This script calculates the Q-Measure as descrived in the DE book. This
% will be the convergence measure divided by the probability of convegence.



function [] = QMeasure (label, results, colFE, colCFV)% Inputs: label --> mark or name to be given to the results so it can be uploaded; results --> results from AMIGO2 simulation to be analysed; colFE --> column in the array were the number of function evaluations is located; colCFV --> column where the cost function value is

%%%%%%%%%%%%%%%%%%%%%%%%%
%  CONVERGENCE MEASURE  %
%%%%%%%%%%%%%%%%%%%%%%%%%

FE = []; % Empty vector that will include all the values of the function evaluations (all the trials)
epsilon = 0; % Still dont know which should be the value of it

for a = 1:length(results) % For loop that goes through the length of the Matrix introduced
    s = results(a,:); % This element extracts each vector (row) of the matrix
    if s(colCFV) <= epsilon % If statement s that only takes into acount cost function values that are smaler or equal to epsilon
        FE = [FE s(colFE)]; % If the if statement is true then the value of the number of function evaluations will be stored in the vector FE
    end 
end

CM = (sum(FE)/length(FE)); % Calculation of the convergence measure as defined

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PROBABILITY OF CONVERGENCE  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

endv = []; % Empty Matrix that will include the last element of each vector of a matrix (under the assumption that the results from AMIGO will be a matrix where the last element of each vector or at least some of them will be the value of the cost function)
nc = 0; % Number of successful trials
nt = length(results); % Total number of trials, which should be 30
epsilon = 0; % Still dont know which should be the value of it

for a = 1:length(results) % For loop that goes through the length of the Matrix introduced
    s = results(a,:); % This element extracts each vector (row) of the matrix
    e = s(end); % This extracts the last element for the vector the loop is on
    endv = [endv e]; % This ads the las element of each vector of the matrix to a new vector that will include this values
end  

for c = 1:length(endv) % For loop that goes through the length of endv
    if endv(c) <= epsilon % If statement s that only takes into acount cost function values that are smaler or equal to epsilon
        nc = nc+1; % If this is true this ads 1 to the value of nc, where at the end will have the total number of succesfull trials
    end
end
    
p = (nc/nt)*100; % Calculation of the probability of convergence or average of converged ones in respect to the total


%%%%%%%%%%%%%%%
%  Q-MEASURE  %
%%%%%%%%%%%%%%%

Qm = CM/p;


mkdir(strcat('Q_Measure',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Makes a new  folder that will contain the results
cd(strcat('Q_Measure',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Changes directory to the new folder

save([label,'-Q_Measure'],'Qm'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Q_Measure is ', num2str(Qm)]); % This just prints the result on the screan to check it

label2 = categorical({['Q-Measure-',label]}); % Labels for the bar charts
bar(label2,Qm ,'r','FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5); % Bar chart that shows in one bar the average number of function evaluations used and in the other bar the maximum boundary assigned to all the algorithm for a visual representation of how many FE an algorithm needed in average in respect to the total

cd ..\; % Back to the original directory



end