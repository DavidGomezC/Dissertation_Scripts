%% From CPU Time to Function Evaluations (DE)

% This function extracts the convergence curve from AMIGO (designed for
% Differential Evolution results), takes the vector with the different
% CPU times and makes a corresponding vector with the values of function
% evaluations at each time. As input takes the name of the directory 
% where the AMIGO2 results are.

function [] = TimeToFunctEvalDE (directory)

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

% For loop that takes all the names of the files from list2 (list with
% desired file names) one by one, loads the file and takes the vector
% correspondent to the CPU time from the convergence curve result and
% stores all of them into the empty vector CPUt. This is because a priory
% all trials will have the same length and including them in an empty
% vector is a check to see if it is true for all (otherwise it will through
% an error. It also takes all the maximum number of function evaluations
% used in every trial and store it into the empty vector mfe.
CPUt = [];
mfe = [];
for x=1:length(list2(:,1))
    load(list2(x,:));
    t = oed_results{1}.nlpsol.conv_curve(:,1);
    fe = oed_results{1}.nlpsol.nfeval;
    CPUt = [CPUt, t];
    mfe = [mfe, fe];
end

% Variable with the length of each vector (same for all) to use as a third
% argument to the function linspace().
l = length(CPUt);
% Average of all the maximum number of function evaluations to be used as a
% second argument for the function linspace().
MFE = sum(mfe)/length(mfe);

% Calculation of the vector with the corrspondent number of function
% evaluations at each time approximately. As an initial point 800 is
% selected since inspection in eSS results indicates that 100 seconds is
% arround 800 function evaluations approximately. 
FE = linspace(800, MFE, l);

% Saves the result of FE using the function input directory
% as a differentiation between all the results generated.
save('FunctionEvaluationsVector', 'FE');

% Back to the original directory
cd ..\Functions;

end