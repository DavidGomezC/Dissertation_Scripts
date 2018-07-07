function [] = AverageOfConvergence (directory, label) % Inputs: label --> mark or name to be given to the results so it can be uploaded; results --> results from AMIGO2 simulation to be analysed

cd ..\;
cd (directory);

SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

endv = []; 
nc = 0; 
nt = length(list2(:,1)); 
epsilon = -1e40; % Provisional for the test

for x=1:length(list2(:,1))
    load(list2(x,:));
    f = oed_results{1}.nlpsol.fbest;
    endv = [endv f];
end

for c = 1:length(endv) % For loop that goes through the length of endv
    if endv(c) <= epsilon % If statement s that only takes into acount cost function values that are smaler or equal to epsilon
        nc = nc+1; % If this is true this ads 1 to the value of nc, where at the end will have the total number of succesfull trials
    end
end

p = (nc/nt)*100; % Calculation of the probability of convergence or average of converged ones in respect to the total

cd ..\; 
cd ('Performance_Measures_New');

cd(strcat('AverageOfConvergence')); % Changes directory to the new folder

save([label,'-ProbabilityOfConvergence'],'p'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Probability of convergence is ', num2str(p),'%']); % This just prints the result on the screan to check it

GrafRep = [p (100-p)]; % Vector with the percentage of converged ones and not converged ones
labels = {['Converged --> ',num2str(p),'%'],['Not Converged --> ', num2str(100-p),'%']}; % Labels for the pie chart
pie3(GrafRep,labels); % Pie chart for a fast visual annalysis of the perfentage that has converged and the one that not

cd ..\

end
