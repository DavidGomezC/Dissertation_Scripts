function [] = MeanFunctEval (directory, label)% Inputs: label --> mark or name to be given to the results so it can be uploaded; results --> results from AMIGO2 simulation to be analysed; colFE --> column in the array were the number of function evaluations is located; maxFE --> maximum value of function evaluations that where conceded to all the simulations; 

maxFV = -5.2159e+42;
FE = []; % Empty vector that will include all the values of the function evaluations (all the trials)

cd ..\;
cd (directory);

SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

cf = {};
fe = {};
for x=1:length(list2(:,1))
    load(list2(x,:));
    a = oed_results{1}.nlpsol.f;
    b = oed_results{1}.nlpsol.neval;
    cf{x} = a;
    fe{x} = b;
end

indx={};
for x=1:length(cf)
    r = cf{x}-maxFV;
    s = abs(floor(log10(maxFV)));
    indx{x} = r < 0.1*10^s;
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
    else
        FE = [FE, 300000];
    end
end



cd ..\; 
cd ('Performance_Measures_New');




MFE = (sum(FE)/length(FE)); % Calculation of the mean number of function evaluations as defined

cd(strcat('MeanFunctionEvaluations')); % Changes directory to the new folder

save([label,'-MeanFuncEval'],'MFE'); % This saves the value of the average number of fuunction evaluations in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Average number of function evaluations is ', num2str(MFE)]); % This just prints the result on the screan to check it

label = categorical({'MFE','Maximum Boundary Assigned'}); % Labels for the bar charts
bar(label,[MFE 300000],'r','FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5); % Bar chart that shows in one bar the average number of function evaluations used and in the other bar the maximum boundary assigned to all the algorithm for a visual representation of how many FE an algorithm needed in average in respect to the total

cd ..\; % Back to the original directory


end