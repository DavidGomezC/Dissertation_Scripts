function [] = Quantiles2 (directory, desiredT)

%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION EVALUATIONS  %
%%%%%%%%%%%%%%%%%%%%%%%%%
maxFV = -5.2159e+34;
FE = [];

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
        FE = [FE, fe{x}(end)];
    end
end

cd ..\; 
cd ('Performance_Measures_New');

med = quantile(FE,0.5); % Takes the median of the dataset

disp(['The median is: ', num2str(med)]); 

p = linspace(0,1, 10000); % Linear space of probabilities from 0 to 1 taking 100 steps (default)
PDT =[]; % Probability(s) to the desired time 
for x=p   % For loop to calculate the probability associated to the quantile with our desired function evaluation value
    if round(quantile(FE,x),2,'significant') == round(desiredT,2,'significant')
        PDT = [PDT, x];
    end
end

PDT = sum(PDT)/length(PDT); % Average of the probabilities

disp(['The probability to converge to a cost function with exponential ', num2str(desiredT), ...
    ' is: ', num2str(PDT)]);

cd(strcat('Quantiles2')); % Changes directory to the new folder

save([directory, '-', num2str(desiredT),'-Quantiles'],'PDT'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.

cd ..\; % Back to the original directory

end