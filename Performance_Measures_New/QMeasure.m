function [] = QMeasure (directory, label)% Inputs: label --> mark or name to be given to the results so it can be uploaded; results --> results from AMIGO2 simulation to be analysed; colFE --> column in the array were the number of function evaluations is located; colCFV --> column where the cost function value is

%%%%%%%%%%%%%%%%%%%%%%%%%
%  CONVERGENCE MEASURE  %
%%%%%%%%%%%%%%%%%%%%%%%%%

FE = []; % Empty vector that will include all the values of the function evaluations (all the trials)
epsilon = -4.3e32; % Provisional for the test
maxFV = -5.2159e+42;
endv = []; % This is for the probability

cd ..\;
cd (directory);

SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end
nc = 0; % This is for the probability
nt = length(list2(:,1)); % This is for the probability
cf = {};
fe = {};
for x=1:length(list2(:,1))
    load(list2(x,:));
    a = oed_results{1}.nlpsol.f;
    b = oed_results{1}.nlpsol.neval;
    cf{x} = a;
    fe{x} = b;
    f = oed_results{1}.nlpsol.fbest; % This is for the probability
    endv = [endv f]; % This is for the probability
end

indx={};
for x=1:length(cf)
    r = cf{x};
    indx{x} = r < epsilon;
end

fir = [];
for x=1:length(indx)
    if find(indx{x}==1)
        i = find(indx{x}==1);
        fir = [fir, i(1)];
    else
        fir = [fir, 0];
    end
end

FE = [];
for x=1:length(fe)
    if fir(x) ~= 0
        FE = [FE, fe{x}(fir(x))];
    end
end

cd ..\; 
cd ('Performance_Measures_New');

CM = (sum(FE)/length(FE)); % Calculation of the convergence measure as defined

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PROBABILITY OF CONVERGENCE  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

cd(strcat('Q_Measure')); % Changes directory to the new folder

save([label,'-Q_Measure'],'Qm'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Q_Measure is ', num2str(Qm)]); % This just prints the result on the screan to check it

cd ..\; % Back to the original directory

end