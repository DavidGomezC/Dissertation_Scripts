% Assesses how long an algorithm is able to keep following the optima since
% its performance might tend to degrade (in terms of fitness quality) as 
% the number of changes becomes larger. This is the slope of the regression
% line of the final accuracy values (recorded at the end of each period). 
% Beta-degradation = quantifies the degradation over time. 
% The loss of solution quality is expressed by a negative slope 
% (bdegradation\0), while a positive value implies the algorithm is able to
% compensate well after changes by producing high quality solutions.
% Nevertheless, since we consider the best value the minimum cost function
% value, a good algorithm will be one with a negative slope since in the
% paper they consider the optima a maximum and for us it is a minimum.
% Aproximation of total accuracy = Beta-degradation * vector of size equal
% to the toal number of periods.
% Two deffinitions of accuracy will be used. The first one is the one used
% in the script Accuracy as it is. The second is how in the papaer they
% define acuracy at each pperiod k in the secction for fitness degradation.
% This is just the average between all the simulations of the best solution
% found so far. So not normalized and the result has to be stored since if
% in a next generation the solution is worst, we still have to keep it as
% the best so far. 

function [] = FitnessDegrad (label, results, worst, best)

%%%%%%%%%%%%%%%%
%  ACCURACY 1  %
%%%%%%%%%%%%%%%%

bsl = best; % Best solution in the llandscape introduced by the user
wks = max(results(:,worst)); % Maximum value of the cost function column (first one)
BS = []; %Empty vector that will include the best solution at each iteration or time

for x=1:length(results)
    bst = min(results(:,x)); % Takes the best value of each time as the minimum of the cost function value for each column
    BS = [BS, bst]; % Puts each best value nto the empty list BS
end

At = []; % Empty vector that will include each acuracy value for each time
for x=1:length(BS)
    a = (Bs(x)-wks)/(bsl-wks); % Claculation of the acuracy value for each time in each loop
    At = [At, a]; % Adition of each valu into At in order
end

%%%%%%%%%%%%%%%%
%  ACCURACY 2  %
%%%%%%%%%%%%%%%%

Rmin = []; % Empty vector that will have the cost function values but just the best ones for each generation
emp=(results*0)+max(results(:)); % Matrix same size of the results but with all its worst value (the maximum)

for x=1:length(results(:,1)) % Nested loop to be able to access in a loop all elements of the matrix (column and row)
    for y = 1:length(results(1,:))
        t = results(x,y); % Value from the results that goes through every value for each row (one by one)
        if t <= min(emp(x,:)) % If the value from the results is smaller than any value of emp (for each row) it takes the minimum into emp
            emp(x,y) = t;
        else                  % If the value from the results is higher then it takes the best value so far
            emp(x,y) =  min(emp(x,:));
        end
    end
end

for x=1:length(emp(1,:))
    v=sum(emp(:,x))/length(emp(:,1));  % Claculation of the average at each generation (which only contains the best of each trial). This is the definition of accuracy they give
    Rmin = [Rmin, v];  % Addition of the results to Rmin
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LINEAR REGRESSION ACCURACY 1  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = [1:length(At)];
y = At;

%r	Regression values for each of the N matrix rows
%m	Slope of regression fit for each of the N matrix rows
%b	Offset of regression fit for each of the N matrix rows

[r,m,b] = regression(x,y);
plotregression(x,y)
title('Linear Regression 1')
xlabel('Periods')
ylabel('Accuracy')

mkdir(strcat('Fitness_Degradation',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Makes a new  folder that will contain the results
cd(strcat('Fitness_Degradation',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Changes directory to the new folder

save([label,'-Beta1'],'m'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Beta degradation 1 is ', num2str(m)]); % This just prints the result on the screan to check it


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR REGRESSION ACCURACY 2  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = [1:length(Rmin)];
y = Rmin;

%r	Regression values for each of the N matrix rows
%m	Slope of regression fit for each of the N matrix rows
%b	Offset of regression fit for each of the N matrix rows

[r2,m2,b2] = regression(x,y);
plotregression(x,y)
title('Linear Regression 2')
xlabel('Periods')
ylabel('Accuracy')

save([label,'-Beta2'],'m2'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Beta degradation 2 is ', num2str(m2)]); % This just prints the result on the screan to check it

cd ..\; % Back to the original directory

end
