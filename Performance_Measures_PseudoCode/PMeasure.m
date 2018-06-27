% Allows us to observe the convergence rate from the point of view of
% variables (dinamics of grouping individuals aroound the optimum). It is 
% supposed that the more densely individuals are located, the better 
% algorithm convergence is intended.Based on the radius of population
% gauging. Euclidian distance between the centre of a population and the 
% individual farthest from it. The center of a population is the average
% vector of all individuals. The problem is that AMIGO2 the value for each
% member of the population at each function evaluation is not stored which
% makes it not possible to use the measure as descrived in the DE book. 

% Hypothesis --> We can modify it as a mesure of convergence rate from the 
% point of view of each stocstic iteration. Population will refear to the number of trials then. 
% So, instead of the average
% vector of individuals we can use the average cost value per function
% evaluation, so instead of per population members we consider per trials.
% The center of the population will be average cost function value for the
% total amount of trials. Hence P-measure will be the euclidean distance
% between the center of a cost value population at each function evaluation
% and the individual fardest from it. 

function [] = PMeasure (label, results)

%for x=1:length(results(:,(length(results)))) % For loop that goes through the length of the colums of the matrix (not vector)

ind=[]; % Empty vector of individuals

for x=1:length(resuts) % For loop that goes trough the length of the vectors in the matrix
    col=results(:,x).'; % This takes each column and makes it a vector
    ind = [ind;col]; % This makes the new matrix which is actualy the transpose of the original one where each vector is one individual (function evaluation)
end

Op=[]; % Empty vector with the center of the population
for x=1:length(ind) % For loop that goes through the length of the vector of individuals
    o=(sum(ind(x,:)))/length(ind(x,:)); % Calculation of the center of the population for each function evaluation
    Op=[Op, o]; % Addition of the different values of the center of the population into one vector
end

dis = []; % Empty vector that will include all the euclidean distances

for x=1:length(Op) % For loop that goes through the length of the vector with the center of the population
    v = []; % Vector that will include the individual distances for each loop that will end up in the matrix dis as different vectors
    for y=1:length(ind(1,:)) % Nested for loop that goes through each row of the matrix of individuals
        r = abs(ind(x,y)-Op(x)); % Calculation of the individual Euclidean distances
        v = [v r]; % Addition of the individual distances for one group into the vector v
    end
    dis = [dis; v]; % Addition of the individual distances to a matrix that will store them all
end

Pm = []; % Empty vector that will include all the P-measure values for each function evaluation or time

for x=1:length(dis) % For loop that goes through the length of the matrix dis with all the distances
    m=max(dis(x,:)); % For each function evaluation it takes the maximum value of the euclidean distance
    Pm = [Pm;m]; % It ads each P-measure value into the empty vector Pm
end

mkdir(strcat('P-Measure',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Makes a new  folder that will contain the results
cd(strcat('P-Measure',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Changes directory to the new folder

save([label,'-PMeasure'],'Pm'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['The P-measure is ', num2str(Pm)]); % This just prints the result on the screan to check it

cd ..\; % Back to the original directory


end
