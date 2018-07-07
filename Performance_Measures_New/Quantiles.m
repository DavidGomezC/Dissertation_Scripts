function [] = Quantiles (directory, desiredT)

%%%%%%%%%%%%%%%%%%%%%%%%
% COST FUNCTION VALUE  %
%%%%%%%%%%%%%%%%%%%%%%%%

cd ..\;
cd (directory);

SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

vect = []; 
for x=1:length(list2(:,1))
    load(list2(x,:));
    f = oed_results{1}.nlpsol.fbest;
    vect = [vect f];
end

cd ('..\Performance_Measures_New');

vect = log10(-vect);

med = quantile(vect,0.5); % Takes the median of the dataset

disp(['The median is: ', num2str(med)]); 

p = linspace(0,1, 10000); % Linear space of probabilities from 0 to 1 taking 100 steps (default)
PDT =[]; % Probability(s) to the desired time 
for x=p   % For loop to calculate the probability associated to the quantile with our desired function evaluation value
    if round(quantile(vect,x)) == round(desiredT)
        PDT = [PDT, x];
    end
end

PDT = sum(PDT)/length(PDT); % Average of the probabilities

disp(['The probability to converge to a cost function with exponential ', num2str(desiredT), ...
    ' is: ', num2str(PDT)]);

cd(strcat('Quantiles')); % Changes directory to the new folder

save([directory, '-', num2str(desiredT),'-Quantiles'],'PDT'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.

cd ..\; % Back to the original directory

end