% This script calculates the average of convergence for a global optimizer
% that has been run X times which is the probability of convergence.

% TEST HYPOTHESIS: introduce my deffinition of probability of convergence
% --> Use the cost function values and average it to the total number
% of trials, then rest it to the best value (converged) to obtain the
% average distance to the best value. The probability will be 1 -
% (Distance/best) where if the cost function stays at 0 then the distance
% will be the same as the best value and hence probability of convergence 0
% while if distance is 0 probability of convergence is 1. This allos to
% consider elements that get close to the optimum but they do not reach it
% (usual tendency of stochastic algorithms). 

function [] = AverageOfConvergence2 (label, results) % Inputs: label --> mark or name to be given to the results so it can be uploaded; results --> results from AMIGO2 simulation to be analysed 

endv = []; % Empty Matrix that will include the last element of each vector of a matrix (under the assumption that the results from AMIGO will be a matrix where the last element of each vector or at least some of them will be the value of the cost function)
BestValue = 1; % Still dont know which should be the value of it (best value obtained by the cost function, the one defined as converged to optimal solution)

for a = 1:length(results) % For loop that goes through the length of the Matrix introduced
    s = results(a,:); % This element extracts each vector (row) of the matrix
    e = s(end); % This extracts the last element for the vector the loop is on
    endv = [endv e]; % This ads the las element of each vector of the matrix to a new vector that will include this values
end

AvOFV = (sum(endv))/length(endv); % Average of the Objective Function values at the last function evaluation

dis = abs(BestValue-AvOFV); % Distance measure from the average to the optimum

p = (1 - abs(dis/BestValue))*100; % Calculation of the probability of a stochastic algorithm to get close to the best value (convergence)

mkdir(strcat('AverageOfConvergence2',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Makes a new  folder that will contain the results
cd(strcat('AverageOfConvergence2',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Changes directory to the new folder

save([label,'-ProbabilityOfConvergence2'],'p'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Probability of convergence is ', num2str(p), '%']); % This just prints the result on the screan to check it


GrafRep = [p (100-p)]; % Vector with the percentage of converged ones and not converged ones
labels = {['Converged --> ',num2str(p),'%'],['Not Converged --> ', num2str(100-p),'%']}; % Labels for the pie chart
pie3(GrafRep,labels); % Pie chart for a fast visual annalysis of the perfentage that has converged and the one that not

cd ..\; % Back to the original directory

end
