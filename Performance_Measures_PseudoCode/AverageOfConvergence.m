% This script calculates the average of convergence for a global optimizer
% that has been run X times which is the probability of convergence
% descrived in the DE book, Succesful trials to total trials. This is
% performed as a first evalutation to see if it is worth to proceed with
% analysis or not (if this probability is 0 there is no sense to progrees).
% Still need to find a way to calculate epsilon as the boundary for
% convergence. 

function [] = AverageOfConvergence (label, results) % Inputs: label --> mark or name to be given to the results so it can be uploaded; results --> results from AMIGO2 simulation to be analysed
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

mkdir(strcat('AverageOfConvergence',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Makes a new  folder that will contain the results
cd(strcat('AverageOfConvergence',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Changes directory to the new folder

save([label,'-ProbabilityOfConvergence'],'p'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Probability of convergence is ', num2str(p),'%']); % This just prints the result on the screan to check it

GrafRep = [p (100-p)]; % Vector with the percentage of converged ones and not converged ones
labels = {['Converged --> ',num2str(p),'%'],['Not Converged --> ', num2str(100-p),'%']}; % Labels for the pie chart
pie3(GrafRep,labels); % Pie chart for a fast visual annalysis of the perfentage that has converged and the one that not

cd ..\; % Back to the original directory

end
