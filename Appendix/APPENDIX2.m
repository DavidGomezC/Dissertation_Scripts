%% APPENDIX 2: Performance Mesures
%
% This appendix contains all the scripts that contain the diverse
% performance measures used in this work.


%% Average Of Convergence 1

% This function calculates the average of convergence for a global 
% optimizer that has been run X times, which is also called probability 
% of convergence. This is the number of succesful trials diveded
% by the number of total trials. As inputs it takes a the name of the
% directory where the AMIGO2 results desired to analyse are and a small
% label that whill be the name of thefile that contains the results of the
% script.

function [] = AverageOfConvergence (directory, label) 

% This goes to the directory where the AMIGO2 results are, in this case in
% a directory that is located in the parent directory where this script is
% located. This can vary deppending on the user.
cd ..\; cd (directory);

% This part of the code just takes the names of all the files of interest
% in the directory so this can be used latter to load all the files.
SN = ls; list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

% Empty vector that will contain the best cost function value of each
% document (results of each trial) loaded
endv = [];
% Number of trials that have converged (starts with 0 and will increase if
% there are any trials that did converged)
nc = 0; 
% Total number of trials
nt = length(list2(:,1));
% Precision range of the convergence
epsilon = -1e40; 

% For loop that takes all the names of the files from list2 (list with
% desired file names) one by one, loads the file and takes the best cost
% function value (the final one for each trial) and stores it in endv.
for x=1:length(list2(:,1))
    load(list2(x,:));
    f = oed_results{1}.nlpsol.fbest;
    endv = [endv f];
end

% For loop that thanks to the if statement compares the cost function value
% to the precission value and if smaller it considers that the tryal has
% converged so increases the variable nc by one.
for c = 1:length(endv) 
    if endv(c) <= epsilon 
        nc = nc+1; 
    end
end

% Calculation of the probability of convergence or average of converged 
% ones in respect to the total
p = (nc/nt)*100; 

% This part of the code goes back to the original directory where this
% script is located. This might vary according to the user.
cd ..\;  cd ('Performance_Measures_New');
% Changes directory to a directory where all the results will be stored
cd(strcat('AverageOfConvergence')); 
% Saves the result of p (probability or average of convergence) using the
% function input label as a differentiation between all the results
% generated.
save([label,'-ProbabilityOfConvergence'],'p'); 
% This just prints the result on the screan to check it
disp(['Probability of convergence is ', num2str(p),'%']); 
% Go back to the directory where this script is located
cd ..\

end

%% Average Of Convergence 2

% This function calculates the average of convergence for a global 
% optimizer that has been run X times, which is also called probability 
% of convergence. As inputs it takes a the name of the
% directory where the AMIGO2 results desired to analyse are and a small
% label that whill be the name of thefile that contains the results of the
% script. As a difference of the Average Of Convergence 1 this script
% calcualtes the probability without the need of a precission value but
% using Euclidean distances. 

function [] = AverageOfConvergence2 (directory, label) 

% Empty vector that will contain the best cost function value of each
% document (results of each trial) loaded
endv = [];
% Best value of the cost function found in total (considering all
% simulations and each of the trials)
BestValue = -5.2159e+42; 

% This goes to the directory where the AMIGO2 results are, in this case in
% a directory that is located in the parent directory where this script is
% located. This can vary deppending on the user.
cd ..\; cd (directory);

% This part of the code just takes the names of all the files of interest
% in the directory so this can be used latter to load all the files.
SN = ls;  list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

% For loop that takes all the names of the files from list2 (list with
% desired file names) one by one, loads the file and takes the best cost
% function value (the final one for each trial) and stores it in endv.
for x=1:length(list2(:,1))
    load(list2(x,:));
    f = oed_results{1}.nlpsol.fbest;
    endv = [endv f];
end

% Average of the Objective Function values at the last function evaluation
AvOFV = (sum(endv))/length(endv); 
% Distance measure from the average to the optimum
dis = abs(BestValue-AvOFV); 
% Calculation of the probability of a stochastic algorithm to get close 
% to the best value (convergence)
p = (1 - abs(dis/BestValue))*100; 

% This part of the code goes back to the original directory where this
% script is located. This might vary according to the user.
cd ..\;  cd ('Performance_Measures_New');
% Changes directory to a directory where all the results will be stored
cd(strcat('AverageOfConvergence2')); 
% Saves the result of p (probability or average of convergence) using the
% function input label as a differentiation between all the results
% generated.
save([label,'-ProbabilityOfConvergence2'],'p'); 
% This just prints the result on the screan to check it
disp(['Probability of convergence is ', num2str(p), '%']);
% Go back to the directory where this script is located
cd ..\; 

end

%% Mean Number of Function Evaluations

% This function calculates the mean (average) number of function
% evaluations that a simulation needed to converge to the global minimum
% with a determined precision. As inputs it takes a the name of the
% directory where the AMIGO2 results desired to analyse are and a small
% label that whill be the name of thefile that contains the results of the
% script.  

function [] = MeanFunctEval (directory, label)

% Best value of the cost function found in total (considering all
% simulations and each of the trials)
maxFV = -5.2159e+42;
% Empty vector that will include all the values of the function 
% evaluations (all the trials)
FE = []; 

% This goes to the directory where the AMIGO2 results are, in this case in
% a directory that is located in the parent directory where this script is
% located. This can vary deppending on the user.
cd ..\; cd (directory);

% This part of the code just takes the names of all the files of interest
% in the directory so this can be used latter to load all the files.
SN = ls;  list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

% Empty cells that will include the cost function values (cf) and function
% evaluation values (fe) for all the trials of one simulation. It has to be
% stored in cells since for each trial the number of cost function values
% (and function evaluations) is different on the results, it only stores a
% value if it is better than the previous one).
cf = {}; fe = {};
% For loop that takes all the names of the files from list2 (list with
% desired file names) one by one, loads the file and takes all the cost
% function and function evaluations to store them in the cells cd, fe.
for x=1:length(list2(:,1))
    load(list2(x,:));
    a = oed_results{1}.nlpsol.f;
    b = oed_results{1}.nlpsol.neval;
    cf{x} = a; fe{x} = b;
end

% This for loop works as the precision value from the script for the 
% average of convergence. It substracts the cost function value to all the
% cost function values of the cell cf and if the results is smaller than
% 0.1 for teen to the power of the same exponential as the best value it
% will sore the position as 1, and if not as 0. This has to be done to know
% the index at which we can consider that the cost function has converged.
indx={};
for x=1:length(cf)
    r = cf{x}-maxFV;
    s = abs(floor(log10(maxFV)));
    indx{x} = r < 0.1*10^s;
end

% Empty vector that will store the index from the original vector for each
% trial at which the cost function converged. 
fir = [];
% For loop that continues the work of the previous one. This just compares
% for each trial if the index is 0 or 1. If in a trial there is an index
% stored as 1 (or more) the number of the first index is stored in fir. If
% there is no index, then the index will be 0.
for x=1:length(indx)
    if find(indx{x}==1)
        i = find(indx{x}==1);
        fir = [fir, i(1)];
    else
        fir = [fir, 0];
    end
end

% Empty vector that will store all the function evaluation values
FE = [];
% For loop that continues the work of the previous one. This takes the
% index stored in fir and extracts the number of function evaluations
% stored in that index for a specific trial. If there was no index (0) then
% the value taken is the maximum number of function evaluations. 
for x=1:length(fe)
    if fir(x) ~= 0
        FE = [FE, fe{x}(fir(x))];
    else
        FE = [FE, 300000];
    end
end

% Calculation of the mean number of function evaluations as defined
MFE = (sum(FE)/length(FE)); 

% This part of the code goes back to the original directory where this
% script is located. This might vary according to the user.
cd ..\;  cd ('Performance_Measures_New');
% Changes directory to a directory where all the results will be stored
cd(strcat('MeanFunctionEvaluations'));
% Saves the result of MFE (Mean function evaluations) using the
% function input label as a differentiation between all the results
% generated.
save([label,'-MeanFuncEval'],'MFE'); 
% This just prints the result on the screan to check it
disp(['Average number of function evaluations is ', num2str(MFE)]); 
% Go back to the directory where this script is located
cd ..\; 

end

%% Quality Measure

% This function calculates the Quality Measure of an algorithm. As inputs 
% it takes a the name of the directory where the AMIGO2 results desired 
% to analyse are, a small label that whill be the name of the file that
% contains the results of the script, the type of desired probability of
% convergence desired to be used and the label that this was stored with.

function [] = QMeasure (directory, label, provability, labelP)

%%%%%%%%%%%%%%%%%%%%%%%%%
%  CONVERGENCE MEASURE  %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Empty vector that will include all the values of the function evaluations (all the trials)
FE = []; 
% Precision range of the convergence
epsilon = -4.3e32;
% Best value of the cost function found in total (considering all
% simulations and each of the trials)
maxFV = -5.2159e+42;

% This goes to the directory where the AMIGO2 results are, in this case in
% a directory that is located in the parent directory where this script is
% located. This can vary deppending on the user.
cd ..\; cd (directory);

% This part of the code just takes the names of all the files of interest
% in the directory so this can be used latter to load all the files.
SN = ls;  list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

% Empty cells that will include the cost function values (cf) and function
% evaluation values (fe) for all the trials of one simulation. It has to be
% stored in cells since for each trial the number of cost function values
% (and function evaluations) is different on the results, it only stores a
% value if it is better than the previous one).
cf = {}; fe = {};
% For loop that takes all the names of the files from list2 (list with
% desired file names) one by one, loads the file and takes all the cost
% function and function evaluations to store them in the cells cd, fe.
for x=1:length(list2(:,1))
    load(list2(x,:));
    a = oed_results{1}.nlpsol.f;
    b = oed_results{1}.nlpsol.neval;
    cf{x} = a; fe{x} = b;
end

% This for loop works as the precision value from the script for the 
% average of convergence. It compares the cost function values to the it
% precision value and if it is smaller it will sore the position as 1, 
% and if not as 0. This has to be done to know the index at which we can 
% consider that the cost function has converged.
indx={};
for x=1:length(cf)
    r = cf{x};
    indx{x} = r < epsilon;
end

% Empty vector that will store the index from the original vector for each
% trial at which the cost function converged. 
fir = [];
% For loop that continues the work of the previous one. This just compares
% for each trial if the index is 0 or 1. If in a trial there is an index
% stored as 1 (or more) the number of the first index is stored in fir. If
% there is no index, then the index will be 0.
for x=1:length(indx)
    if find(indx{x}==1)
        i = find(indx{x}==1);
        fir = [fir, i(1)];
    else
        fir = [fir, 0];
    end
end

% Empty vector that will store all the function evaluation values
FE = [];
% For loop that continues the work of the previous one. This takes the
% index stored in fir and extracts the number of function evaluations
% stored in that index for a specific trial. If there was no index (0) then
% the value taken is the maximum number of function evaluations. 
for x=1:length(fe)
    if fir(x) ~= 0
        FE = [FE, fe{x}(fir(x))];
    end
end

% This part of the code goes back to the original directory where this
% script is located. This might vary according to the user.
cd ..\;  cd ('Performance_Measures_New');

% Calculation of the convergence measure as defined
CM = (sum(FE)/length(FE)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PROBABILITY OF CONVERGENCE  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If the user wants to ue the probability calculated by the script
% AverageOfConvergence then it will move to the directory and load the
% desired result. It is the same if the AverageOfConvergence2 is the one
% desired to be used
if provability == 1
    cd ('AverageOfConvergence');
    load([labelP, '-ProbabilityOfConvergence']);
elseif provability == 2
    cd ('AverageOfConvergence2');
    load([labelP, '-ProbabilityOfConvergence2']);
end

% Go back to the directory where this script is located
cd ..\;

%%%%%%%%%%%%%%%
%  Q-MEASURE  %
%%%%%%%%%%%%%%%

% Calculation of the Q-measure as defined
Qm = CM/p;

% Changes directory to a directory where all the results will be stored
cd(strcat('Q_Measure')); 
% Saves the result of Qm (Quality Measure) using the
% function input label as a differentiation between all the results
% generated.
save([label,'-Q_Measure'],'Qm');  
% This just prints the result on the screan to check it
disp(['Q_Measure is ', num2str(Qm)]);
% Go back to the directory where this script is located
cd ..\; % Back to the original directory

end

%% Robustness Measure

% This function calculates the Robustness Measure of an algorithm as 
% defined in this work, incuding all the modifications. As inputs 
% it takes a the name of the directory where the desired measure to be 
% used is located results desired to analyse are and a small label that 
% whill be the name of the file that contains the results of the script.

function [] = RMeasure (measure, label)

%%%%%%%%%%%%%%
%  DISTANCE  %
%%%%%%%%%%%%%%

% Changes the directory where theperformance measures to analyse are
% located. The structure will vary according to the user.
cd (measure);

% This part of the code just takes the names of all the files of interest
% in the directory so this can be used latter to load all the files.
SN = ls; list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),label)) 
        list2=[list2; SN(x,:)]; 
    end
end

% Empty vectors that will store the performance measures of all the trials
% desired to use (vper) and extract the best value according to the measure
% used.
vper = []; best = [];
% For loop that thanks to the if statement loads and stores all the
% performance measure values into vper and extracts the best one of them
% into best.
for x=1:length(list2(:,1))
    load(list2(x,:));
    if strcmp(measure, 'Q_Measure')
        vper = [vper, Qm];
        best = min(vper);
    elseif strcmp(measure, 'AverageOfConvergence')
        vper = [vper, p];
        best = max(vper);    
    elseif strcmp(measure, 'MeanFunctionEvaluations')
        vper = [var, MFE];
        best = min(vper);
    end
end
% Go back to the original directory.
cd ..\;
% Mean number of the results from the performance measure used to 
% calculate the robustness
mn = sum(vper)/length(vper); 
% Euclidean distance from the best value to the mean of values
dis=abs(best-mn); 
% Printing of the results on the screan
disp(['Best value is: ',num2str(best)]); 
disp(['The mean value is: ',num2str(mn)]);
disp(['The difference between the two values is: ',num2str(dis)]);
fprintf('\n*\nNEXT\n*\n');
 
%%%%%%%%%%%%%%%%%%%%%%%%
%  STANDARD DEVIATION  %
%%%%%%%%%%%%%%%%%%%%%%%%

% Empty vector that will store part of the standard deviation formula to be
% summed.
par=[]; 
% For loop that goes through all the results in vper and calculates the
% square of the rest between each value and the mean which will be stored 
% in par and will have to be summed and averaged.
for c=1:length(vper) 
    e=(vper(c)-mn)^2; 
    par = [par e]; 
end
% Calculation of the standard deviation
sd = sqrt(sum(par)/length(par)); 

% Changes directory to a directory where all the results will be stored
cd(strcat('R-Measure')); 
% Saves the result of sd (standard deviation) using the function input 
% label as a differentiation between all the results generated.
save([label,'-Standard Deviation'],'sd'); 
% This just prints the result on the screan to check it
disp(['The standard deviation is ', num2str(sd)]); 
fprintf('\n*\nNEXT\n*\n');

%%%%%%%%%%%%%%
%  VARIANCE  %
%%%%%%%%%%%%%%

% Calculation of the Variance
vr = (sd)^2; 

% Saves the result of vr (variance) using the function input 
% label as a differentiation between all the results generated.
save([label,'-Variance'],'vr'); 
% This just prints the result on the screan to check it
disp(['The variance is ', num2str(vr)]); 
fprintf('\n*\nNEXT\n*\n');

%%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFIED DEVIATION  %
%%%%%%%%%%%%%%%%%%%%%%%%

% All this part of code is exactly the same as the part regarding standard
% deviation but instead of using the mean value it uses the best value.
parm=[]; 
for c=1:length(vper) 
    e=(vper(c)-best)^2; 
    parm = [parm e]; 
end
msd = sqrt(sum(parm)/length(parm)); 

save([label,'-ModifiedStandardDeviation'],'msd'); 
disp(['The modified standard deviation is ', num2str(msd)]); 
fprintf('\n*\nNEXT\n*\n');

%%%%%%%%%%%%%%%%%%%%%%%
%  MODIFIED VARIANCE  %
%%%%%%%%%%%%%%%%%%%%%%%

% All this part of code is exactly the same as the part regarding variance
% but instead of using the mean value it uses the best value.
mvr = (msd)^2; 

save([label,'-ModifiedVariance'],'mvr');
disp(['The modified variiancen is ', num2str(mvr)]); 
% Back to the original directory
cd ..\; 
end

%% Population Measure

% This function calculates Population Measure of an algorithm. As inputs 
% it takes a the name of the directory where the AMIGO2 results desired 
% to analyse are and a small label that whill be the name of the file that
% contains the results of the script.

function [] = PMeasure (directory, label)

% Empty cell of individuals
ind={}; 

% This goes to the directory where the AMIGO2 results are, in this case in
% a directory that is located in the parent directory where this script is
% located. This can vary deppending on the user.
cd ..\; cd (directory);

% This part of the code just takes the names of all the files of interest
% in the directory so this can be used latter to load all the files.
SN = ls;  list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

% Empty cells that will include the cost function values (cf) and function
% evaluation values (fe) for all the trials of one simulation. It has to be
% stored in cells since for each trial the number of cost function values
% (and function evaluations) is different on the results, it only stores a
% value if it is better than the previous one).
cf = {}; fe = {};
% For loop that takes all the names of the files from list2 (list with
% desired file names) one by one, loads the file and takes all the cost
% function and function evaluations to store them in the cells cd, fe.
for x=1:length(list2(:,1))
    load(list2(x,:));
    s = oed_results{1}.nlpsol.neval;
    c = oed_results{1}.nlpsol.f;
    cf{x} = c; fe{x} = s;
end

% This part of the code goes back to the original directory where this
% script is located. This might vary according to the user.
cd ..\; cd ('Performance_Measures_New');

% x is the length of the x axis (more than 300000 to avoid any errors or
% miss values) and y an empty cell that will contain the cost function
% values
x = 1:300800; y = {};

% This parfor loop takes all the cost function values for each trial and
% generates different vectors where if for an x value there is a cost
% function value, it will include that in the correct index, if not it will
% include a NaN in the index. This has to be done for correct comparison of
% all the trials, since the cost function value obtained in AMIGO is not at
% the same point every time but at the points where a better value is
% obtained. 
parfor q=1:length(cf)
    ind = [];
    for w=1:length(x)
        if ismember(x(w), fe{q})
            ind = find(fe{q}==x(w));
            y{q}(w) = cf{q}(ind);
        else
            y{q}(w) = nan;
        end
    end    
end

% Empty vector that will include all the cost function vectors with all the
% values.
F = [];

% This for loop substitutes each NaN value in the vector for the previous
% value (we assume that the values in between are worst, so we just
% consider the best one at each generation)
for q=1:length(y)
    F = [F; fillmissing(y{q},'previous')]; 
end
% Substitution of remaining NaN to 0 at the begining since there is no
% previous value.
F(isnan(F))=0; 

% Calculation of the center of the population at each function evaluation
OP = sum(F)/length(F(:,1));

% Calculation of the distance to the center of the population at each
% function evaluation for all the trials
dis = abs(F - OP); 

% Selection of the farthest individual to thecenter of the population at
% each function evaluation. Vector is multiplied by -1 since the absolute
% values have been used and for better comparison with the convergence
% curve results average (OP). 
Pm = max(dis)*(-1); 

% Changes directory to a directory where all the results will be stored
cd(strcat('P-Measure')); 
% Saves the result of Pm (Population Measure) using the
% function input label as a differentiation between all the results
% generated.
save([label,'-PMeasure'],'Pm'); 

%%%%% Save Results CSV %%%%%
% This is done to be able to plot the results using python
PMcsv = [1, Pm].'; OPcsv = [2, OP].';
abc = [PMcsv(:), OPcsv(:)]; csvwrite([label,'PMeasure.csv'],abc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Back to the original directory
cd ..\; 

end

%% Accuracy (relative error)

% This function calculates the accuracy of an algorithm. As inputs 
% it takes a the name of the directory where the AMIGO2 results desired 
% to analyse are which will also be the name of the result files generated.

function [] = Accuracy (directory) 

% ACCESS BEST SOLUTION
cd('..\Functions\Bestf'); load('BestTotal_f_Global'); bsl = m;

% ACCESS WORST SOLUTION
cd ('..\Worstf'); load('WorstTotal_f_Global'); wks = w;

% This goes to the directory where the AMIGO2 results are, in this case in
% a directory that is located in the parent directory where this script is
% located. This can vary deppending on the user.
cd ('..\..\'); cd (directory);

% This part of the code just takes the names of all the files of interest
% in the directory so this can be used latter to load all the files.
SN = ls;  list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

% Empty cells that will include the cost function values (cf) and function
% evaluation values (fe) for all the trials of one simulation. It has to be
% stored in cells since for each trial the number of cost function values
% (and function evaluations) is different on the results, it only stores a
% value if it is better than the previous one).
cf = {}; fe={};
% For loop that takes all the names of the files from list2 (list with
% desired file names) one by one, loads the file and takes all the cost
% function and function evaluations to store them in the cells cd, fe.
for x=1:length(list2(:,1))
    load(list2(x,:));
    f = oed_results{1}.nlpsol.f;
    s = oed_results{1}.nlpsol.neval;
    cf{x} = f; fe{x} = s;
end

% This part of the code goes back to the original directory where this
% script is located. This might vary according to the user.
cd ..\;  cd ('Performance_Measures_New');

% x is the length of the x axis (more than 300000 to avoid any errors or
% miss values) and y an empty cell that will contain the cost function
% values
x = 1:300800; cft = {};

% This parfor loop takes all the cost function values for each trial and
% generates different vectors where if for an x value there is a cost
% function value, it will include that in the correct index, if not it will
% include a NaN in the index. This has to be done for correct comparison of
% all the trials, since the cost function value obtained in AMIGO is not at
% the same point every time but at the points where a better value is
% obtained. 
parfor q=1:length(cf)
    ind = [];
    for w=1:length(x)
        if ismember(x(w), fe{q})
            ind = find(fe{q}==x(w));
            cft{q}(w) = cf{q}(ind);
        else
            cft{q}(w) = nan;
        end
    end    
end

% Empty vector that will include all the cost function vectors with all the
% values.
F = [];
% This for loop substitutes each NaN value in the vector for the previous
% value (we assume that the values in between are worst, so we just
% consider the best one at each generation)
for q=1:length(cft)
    F = [F; fillmissing(cft{q},'previous')]; 
end
% Substitution of remaining NaN to 0 at the begining since there is no
% previous value.
F(isnan(F))=0; 

% Best solution at each function evaluation
BS = min(F);

% Calculation of accuracy as deffined
At = (BS-wks)/(bsl-wks);
% Changes directory to a directory where all the results will be stored
cd(strcat('Accuracy')); 

%%%%% Save Results CSV %%%%%
xcsv = [0, x].'; ycsv = [1, y].';
acc = [xcsv(:), ycsv(:)]; csvwrite([directory,'Accuracy.csv'],acc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saves the result of At (Accuracy) using the function input directory
% as a differentiation between all the results generated.
save([directory,'-Accuracy'],'At'); 
% This just prints the result on the screan to check it
disp(['Acuracy al the last time is time is ', num2str(At(end))]); 
% Back to the original directory
cd ..\; 

end

%% Stability

% This function calculates the stability of an algorithm using its 
% accuracy as reference. As inputs it takes a the name of the directory 
% where the Accuracy results are stored.

function [] = Stability (directory)

%%%%%%%%%%%%%%
%  ACCURACY  %
%%%%%%%%%%%%%%

% Access and load the accuracy results generated with the script Accuracy.
% Its location will vary according to where the user has located them. 
Accuracy(directory); cd ('Accuracy\');
load([directory,'-Accuracy']); cd ..\;

%%%%%%%%%%%%%%%
%  STABILITY  %
%%%%%%%%%%%%%%%

% Empty vector that will contain all the values for the stability over generations
Sg = []; 

% For loop that calculates the stability value at each function evaluation
% where all the negative values will be considered as 0 and saved in the
% vector Sg.
for x=1:length(At)
    % To avoid error with the first value 
    if x-1 ~= 0 
        s = At(x)-At(x-1); 
        if s < 0 
            s = 0;
        end
        Sg = [Sg, s]; 
    end
end

% Changes directory to a directory where all the results will be stored
cd(strcat('Stability')); 

%%%%% Save Results CSV %%%%%
st = [1, Sg].'; csvwrite([directory,'Stability.csv'],st);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saves the result of Sg (Stability) using the function input directory
% as a differentiation between all the results generated.
save([directory,'-Stability'],'Sg'); 

cd ..\; % Back to the original directory

end

%% Fitness Degradation

% This function calculates the fitness degradation of an algorithm 
% using its accuracy as reference. As inputs it takes a the name 
% of the directory where the results will be stored.

function [] = FitnessDegrad (directory)

%%%%%%%%%%%%%%%
%  ACCURACY   %
%%%%%%%%%%%%%%%

% Access and load the accuracy results generated with the script Accuracy.
% Its location will vary according to where the user has located them. 
Accuracy(directory); cd ('Accuracy\');
load([directory,'-Accuracy']); cd ..\;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR REGRESSION ACCURACY   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Definition of the x and y axis of the scatter plot
x = [1:length(At)]; y = At;
% Calculation of the regression line elements (slope and starting point)
coef = polyfit(x,y,1);
% Function that represents the linear regression
rl1 = coef(1)*x + coef(2);
% Slope
beta = coef(1);

% Saves the result of beta (slope for At) using the function input 
% directory as a differentiation between all the results generated.
save([directory,'-Beta'],'beta'); 
% This just prints the result on the screan to check it
disp(['Beta degradation is ', num2str(beta)]); 

%%%%% Save Results CSV %%%%%
A1 = [1, y].'; R1 = [2, rl1].';
FD = [A1(:), R1(:)];
csvwrite([directory,'FitnessDegradation.csv'],FD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Back to the original directory
cd ..\; 

end

%% Area Between Curves

% This function calculates the area between two accuracy curves. As 
% inputs it takes a the two accuracy curve results that will be analised. 

function [] = AreaBetweenCurves (resultsA1, resultsA2) 

% Access and load the accuracy results generated with the script Accuracy.
% Its location will vary according to where the user has located them. 
cd('Accuracy');

% This part of the code just takes the names of all the files of interest
% in the directory so this can be used latter to load all the files.
SN = ls;  list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:), resultsA1)) || (contains(SN(x,:),resultsA2))
        list2=[list2; SN(x,:)]; 
    end
end

% For loop that takes all the names of the files from list2 (list with
% desired file names) one by one, loads the file and takes all the two
% accuracy vectors to be analysed. 
curvs = {};
for x = 1:length(list2(:,1))
    load(list2(x,:));
    curvs{x} = At;
end

% Change back to the original directory
cd ..\;
% Define the x axis of the function and plot
x1 = 1:length(curvs{1});  x2 = 1:length(curvs{2});
% Definition of the two curves as y axis 1 and 2
y1 = curvs{1}; y2 = curvs{2};

% Calculation of the area under each curve
pA1 = trapz(x1,y1);  pA2 = trapz(x2,y2);

% Calculation of the area between curves
ABC = pA1 - pA2; 

% Changes directory to a directory where all the results will be stored
cd(strcat('AreaBetweenCurves')); 

% Saves the result of ABC using the function input directory
% as a differentiation between all the results generated.
save([resultsA1,'-', resultsA2,'-AreaBetweenCurves'],'ABC'); 
% This just prints the result on the screan to check it
disp(['The area between curves is ', num2str(ABC)]); 

%%%%% Save Results CSV %%%%%
y1csv = [1, y1].'; y2csv = [2, y2].';
abc = [y1csv(:), y2csv(:)];
csvwrite([resultsA1, resultsA2,'ABC.csv'],abc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Back to the original directory
cd ..\; 

end

%% Quantiles 1

% This function calculates the probability of obtaining at least a
% determined cost function value (desiredT). The second input is the name
% of the directory where the different AMIGO2 results are. 

function [] = Quantiles (directory, desiredT)

%%%%%%%%%%%%%%%%%%%%%%%%
% COST FUNCTION VALUE  %
%%%%%%%%%%%%%%%%%%%%%%%%

% This goes to the directory where the AMIGO2 results are, in this case in
% a directory that is located in the parent directory where this script is
% located. This can vary deppending on the user.
cd ..\; cd (directory);
% This part of the code just takes the names of all the files of interest
% in the directory so this can be used latter to load all the files.
SN = ls;  list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

% For loop that takes all the names of the files from list2 (list with
% desired file names) one by one, loads the file and takes the best cost
% function value (the final one for each trial) and stores it in vect.
vect = []; 
for x=1:length(list2(:,1))
    load(list2(x,:));
    f = oed_results{1}.nlpsol.fbest;
    vect = [vect f];
end

% This part of the code goes back to the original directory where this
% script is located. This might vary according to the user.
cd ('..\Performance_Measures_New');
% Extaction of the exponential for the cost function values (complete cost
% function values generates some problems). 
vect = log10(-vect);
% Takes the median of the dataset
med = quantile(vect,0.5); 
disp(['The median is: ', num2str(med)]); 

% Linear space of probabilities from 0 to 1 taking 10000 steps
p = linspace(0,1, 10000); 

% Location of the value introduced in one of the quantiles and storage of
% the probability associated to it.
PDT =[]; 
for x=p   
    if round(quantile(vect,x)) == round(desiredT)
        PDT = [PDT, x];
    end
end

% Average of the probabilities (this is because we generate a big range of
% quantiles)
PDT = sum(PDT)/length(PDT); 

% Display of the result
disp(['The probability to converge to a cost function with exponential ', num2str(desiredT), ...
    ' is: ', num2str(PDT)]);

% Changes directory to a directory where all the results will be stored
cd(strcat('Quantiles')); 
% Saves the result of PDT (probability) using the function input directory
% as a differentiation between all the results generated.
save([directory, '-', num2str(desiredT),'-Quantiles'],'PDT'); 
% Back to the original directory
cd ..\; 

end

%% Quantiles 2

% This function calculates the probability of obtaining at least the best
% cost function value at a determined number of function evaluations 
% (desiredT). The second input is the name of the directory where the 
% different AMIGO2 results are. 

function [] = Quantiles2 (directory, desiredT)

%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION EVALUATIONS  %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Best cost function value and vector that will store the function
% evaluations
maxFV = -5.2159e+34; FE = [];

% This goes to the directory where the AMIGO2 results are, in this case in
% a directory that is located in the parent directory where this script is
% located. This can vary deppending on the user.
cd ..\; cd (directory);
% This part of the code just takes the names of all the files of interest
% in the directory so this can be used latter to load all the files.
SN = ls;  list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

% Empty cells that will include the cost function values (cf) and function
% evaluation values (fe) for all the trials of one simulation. It has to be
% stored in cells since for each trial the number of cost function values
% (and function evaluations) is different on the results, it only stores a
% value if it is better than the previous one).
cf = {}; fe={};
% For loop that takes all the names of the files from list2 (list with
% desired file names) one by one, loads the file and takes all the cost
% function and function evaluations to store them in the cells cd, fe.
for x=1:length(list2(:,1))
    load(list2(x,:));
    a = oed_results{1}.nlpsol.f;
    b = oed_results{1}.nlpsol.neval;
    cf{x} = a; fe{x} = b;
end

% This for loop works as the precision value from the script for the 
% average of convergence. It substracts the cost function value to all the
% cost function values of the cell cf and if the results is smaller than
% 0.1 for teen to the power of the same exponential as the best value it
% will sore the position as 1, and if not as 0. This has to be done to know
% the index at which we can consider that the cost function has converged.
indx={};
for x=1:length(cf)
    r = cf{x}-maxFV;
    s = abs(floor(log10(maxFV)));
    indx{x} = r < 0.1*10^s;
end

% Empty vector that will store the index from the original vector for each
% trial at which the cost function converged. 
fir = [];
% For loop that continues the work of the previous one. This just compares
% for each trial if the index is 0 or 1. If in a trial there is an index
% stored as 1 (or more) the number of the first index is stored in fir. If
% there is no index, then the index will be 0.
for x=1:length(indx)
    if find(indx{x}==1)
        i = find(indx{x}==1);
        fir = [fir, i(1)];
    else
        fir = [fir, 0];
    end
end

% For loop that continues the work of the previous one. This takes the
% index stored in fir and extracts the number of function evaluations
% stored in that index for a specific trial. If there was no index (0) then
% the value taken is the maximum number of function evaluations. 
for x=1:length(fe)
    if fir(x) ~= 0
        FE = [FE, fe{x}(fir(x))];
    else
        FE = [FE, fe{x}(end)];
    end
end

% This part of the code goes back to the original directory where this
% script is located. This might vary according to the user.
cd ..\;  cd ('Performance_Measures_New');

% Takes the median of the dataset
med = quantile(FE,0.5); 
% Display of the results
disp(['The median is: ', num2str(med)]); 

% Linear space of probabilities from 0 to 1 taking 10000 steps
p = linspace(0,1, 10000); 

% Location of the value introduced in one of the quantiles and storage of
% the probability associated to it.
PDT =[]; 
for x=p   
    if round(quantile(FE,x),2,'significant') == round(desiredT,2,'significant')
        PDT = [PDT, x];
    end
end

% Average of the probabilities (this is because we generate a big range of
% quantiles)
PDT = sum(PDT)/length(PDT); 

% Display of the result
disp(['The probability to converge to a cost function with exponential ', num2str(desiredT), ...
    ' is: ', num2str(PDT)]);

% Changes directory to a directory where all the results will be stored
cd(strcat('Quantiles2')); 
% Saves the result of PDT (probability) using the function input directory
% as a differentiation between all the results generated.
save([directory, '-', num2str(desiredT),'-Quantiles'],'PDT'); 
% Back to the original directory
cd ..\; 

end