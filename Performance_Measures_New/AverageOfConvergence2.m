function [] = AverageOfConvergence2 (directory, label) % Inputs: label --> mark or name to be given to the results so it can be uploaded; results --> results from AMIGO2 simulation to be analysed 

endv = []; % Empty Matrix that will include the last element of each vector of a matrix (under the assumption that the results from AMIGO will be a matrix where the last element of each vector or at least some of them will be the value of the cost function)
BestValue = -5.2159e+42; % Still dont know which should be the value of it (best value obtained by the cost function, the one defined as converged to optimal solution)

cd ..\;
cd (directory);

SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end
for x=1:length(list2(:,1))
    load(list2(x,:));
    f = oed_results{1}.nlpsol.fbest;
    endv = [endv f];
end

AvOFV = (sum(endv))/length(endv); % Average of the Objective Function values at the last function evaluation

dis = abs(BestValue-AvOFV); % Distance measure from the average to the optimum

p = (1 - abs(dis/BestValue))*100; % Calculation of the probability of a stochastic algorithm to get close to the best value (convergence)

cd ..\; 
cd ('Performance_Measures_New');

cd(strcat('AverageOfConvergence2')); % Changes directory to the new folder

save([label,'-ProbabilityOfConvergence2'],'p'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Probability of convergence is ', num2str(p), '%']); % This just prints the result on the screan to check it


GrafRep = [p (100-p)]; % Vector with the percentage of converged ones and not converged ones
labels = {['Converged --> ',num2str(p),'%'],['Not Converged --> ', num2str(100-p),'%']}; % Labels for the pie chart
pie3(GrafRep,labels); % Pie chart for a fast visual annalysis of the perfentage that has converged and the one that not

cd ..\; % Back to the original directory

end
