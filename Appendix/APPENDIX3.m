%% APPENDIX 3: Other Functions
%
% This appendix contains scripts of different functions used in 
% this work to ease different tasks but that are not fully related to 
% the different simulations or performance measures.

%% Best Cost Function Value

% This function extracts the best cost function value for each trial of 
% the different simulations, hence the best local value. As input
% takes the name of the directory where the AMIGO2 results are.

function [] = Bestf (directory)

% This goes to the directory where the AMIGO2 results are, in this case in
% a directory that is located in the parent directory where this script is
% located. This can vary depending on the user.
cd ..\; cd (directory);

% This part of the code just takes the names of all the files of interest
% in the directory so this can be used later to load all the files.
SN = ls; list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

% For loop that takes all the names of the files from list2 (list with
% desired file names) one by one, loads the file and takes the best cost
% function value (the final one for each trial) and stores it in f_val.
f_val = [];
for x=1:length(list2(:,1))
    load(list2(x,:));
    f = oed_results{1}.nlpsol.fbest;
    f_val = [f_val, f];
end

% Obtains the best cost function value (the smallest one).
m = min(f_val);

% This changes the directory to where the results will be stored. This will
% vary depending on the user.
cd ..\; cd ('Functions'); cd ('Bestf');

% Saves the result of m using the function input directory
% as a differentiation between all the results generated.
save(['Best_f_',directory], 'm');
% Back to the original directory
cd ..\

end

%% Best Total Cost Function Value

% This function extracts the best cost function value for all the
% simulations and subsequent trials, hence the best global value. 

function [] = Bestf_total ()

% This changes the directory where the AMIGO2 results are located. This 
% will vary depending on the user. 
cd ..\..\; 
% This part of the code just takes the names of all the directories of 
% interest to use as each one as an input for the function Bestf. This will
% be tagged with the string simul
SN = ls; list1 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'simul')) 
        list1=[list1; SN(x,:)]; 
    end
end

% Move directory to where the function Bestf is located. This will vary
% depending on the user.
cd ('Functions');
% Call of the function Bestf for each of the simulations going through the
% vector list1 which contains the name of all the directories.
for x=1:length(list1(:,1))
    Bestf(list1(x,:));
end
% This changes the directory where the best cost function value for each
% simulation is located. This will vary depending on the user. 
cd ('Bestf');
% This part of the code just takes the names of all the files of interest
% in the directory so this can be used later to load all the files.
SN = ls; list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Best_f')) 
        list2=[list2; SN(x,:)]; 
    end
end

% Load of all the Bestf results sotred into the empty vector Bf
Bf = [];
for x=1:length(list2(:,1))
    load(list2(x,:));
    Bf = [Bf,m];
end

% Obtains the best cost function value (the smallest one).
BF = min(Bf);
% Saves the result of BF using the function input directory
% as a differentiation between all the results generated.
save('BestTotal_f_Global', 'BF');
% Displays the results to check
disp(['The best cost function value (mnimum) is: ',num2str(BF)]);

end

%% Worst Cost Function Value

% This function extracts the worst cost function value for each trial of 
% the different simulations, hence the worst local value. As input
% takes the name of the directory where the AMIGO2 results are.

function [] = Worstf (directory)

% This goes to the directory where the AMIGO2 results are, in this case in
% a directory that is located in the parent directory where this script is
% located. This can vary depending on the user.
cd ..\; cd (directory);

% This part of the code just takes the names of all the files of interest
% in the directory so this can be used later to load all the files.
SN = ls; list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

% For loop that takes all the names of the files from list2 (list with
% desired file names) one by one, loads the file and takes the worst cost
% function value (the first one for each trial) and stores it in f_val.
f_val = [];
for x=1:length(list2(:,1))
    load(list2(x,:));
    f = oed_results{1}.nlpsol.f(1);
    f_val = [f_val, f];
end

% Obtains the worst cost function value (the largest one).
w = max(f_val);

% This changes the directory to where the results will be stored. This will
% vary depending on the user.
cd ..\; cd ('Functions'); cd ('Worstf');
% Saves the result of w using the function input directory
% as a differentiation between all the results generated.
save(['Worst_f_',directory], 'w');
% Back to the original directory
cd ..\

end

%% Worst Total Cost Function Value

% This function extracts the worst cost function value for all the
% simulations and subsequent trials, hence the worst global value. 

function [] = Worstf_total ()

% This changes the directory where the AMIGO2 results are located. This 
% will vary depending on the user. 
cd ..\..\;
% This part of the code just takes the names of all the directories of 
% interest to use as each one as an input for the function Worstf. 
% This will be tagged with the string simul
SN = ls; list1 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'simul')) 
        list1=[list1; SN(x,:)]; 
    end
end

% Move directory to where the function Worstf is located. This will vary
% depending on the user.
cd ('Functions');
% Call of the function Worstf for each of the simulations going through 
% the vector list1 which contains the name of all the directories.
for x=1:length(list1(:,1))
    Worstf(list1(x,:));
end
% This changes the directory where the best cost function value for each
% simulation is located. This will vary depending on the user. 
cd ('Worstf');
% This part of the code just takes the names of all the files of interest
% in the directory so this can be used later to load all the files.
SN = ls;  list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Worst_f')) 
        list2=[list2; SN(x,:)]; 
    end
end

% Load of all the Worstf results sotred into the empty vector Wf
Wf = [];
for x=1:length(list2(:,1))
    load(list2(x,:));
    Wf = [Wf,w];
end

% Obtains the worst cost function value (the largest one).
WF = max(Wf);
% Saves the result of WF using the function input directory
% as a differentiation between all the results generated.
save('WorstTotal_f_Global', 'w');
% Displays the results to check
disp(['The Worst cost function value (mnimum) is: ',num2str(WF)]);

end

%% Convergence Curves

% This script extracts all the convergence curves of a simulation and save
% them in a CSV format to be plotted in Python. As arguments it takes directory,
% which is the name of the directory where the results of AMIGO2 are and
% also a label name for all the resultant CSV files and the directory where
% these will be stored.

function [] = ConvCurvToCSV (label, directory)

% This changes the directory where the AMIGO2 results are located. This 
% will vary depending on the user. 
cd ..\; cd(directory); 

% This part of the code just takes the names of all the directories of 
% interest to use as each one as an input for the function Worstf. 
% This will be tagged with the string simul
SN = ls; list2 =[]; 
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

% Generates a new directory where the resultant CSV files will be stored
mkdir([label, '_CC_CSV']); cd([label, '_CC_CSV']);

% Save the different convergence curves in separated CSV files. This has to
% be done in separated files since the curves will have different vector
% length.
for x=1:length(cf)
    cc = []; cc = [fe{x}; cf{x}];
    csvwrite([label, num2str(x),'.csv'],cc);
end
% Change directory to the original one. This will vary depending on the
% user.
cd ..\..\; cd('Functions');

end

%% Input Profile

% This script extracts all the best input profile of a simulation and 
% saves it in a CSV format to be plotted in Python. As arguments it takes 
% directory, which is the name of the directory where the results 
% of AMIGO2.

function [] = InputProf (directory)

% This changes the directory where the AMIGO2 results are located. This 
% will vary depending on the user. 
cd ..\; cd(directory);
% Sets the best convergence curve value
best = -5e+42;
% This part of the code just takes the names of all the directories of 
% interest to use as each one as an input for the function Worstf. 
% This will be tagged with the string simul
SN = ls;  list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

% For loop that takes all the names of the files from list2 (list with
% desired file names) one by one, loads the file and extracts the best
% input profile for the trial with the best convergence value achieved in
% the simulation set be the variable best.
for x=1:length(list2(:,1))
    load(list2(x,:));
    if oed_results{1}.nlpsol.fbest <= best
        IP = oed_results{1}.nlpsol.vbest;
    end
end

% Repetition of each value 200 times to be plotted accordingly to time and
% stored in the vector u.
I = [0, IP]; u = [];
for x=1:length(IP)
    r = repelem(IP(x),200); u = [u, r];
end

% Generates a new directory where the resultant CSV files will be stored
mkdir('InputProfile'); cd('InputProfile');
% Write the resultant input profile (single values and repeated 200 times)
% in a CSV file
csvwrite('InputP.csv',u.');
csvwrite('InputSing.csv',I.');
% Change directory to the original one. This will vary depending on the
% user.
cd ..\..\; cd('Functions');

end