% In stochastic algoorithms the distribution of solution quality is often
% assymetrical. How this works is that  the set of solution quality or
% results are ordered from smallesst to highest where the median will be
% the quantile with probability 0.5 and all quantiles will have a
% probability p associated. Hence, if a quantile of probability 0.7 has a
% value of 36, then this means that there is AT LEAST a 0.7 probability of
% getting a solution with the quality equal to 36 or better.
% FROM PAPER --> A big advantage of using quantiles is their applicability to
% multiple executions of the same algorithm for a particular
% problem instance. In the common approach, when
% computational resources are restricted and the algorithm is
% independently executed in a sequential or parallel manner, the
% best achieved solution is used as the final result. If the
% algorithm with restricted computational resources achieves
% quantile Qp of solution quality, then for n independent
% executions there is probability 1 – (1 – p)n of getting at least
% once the solution of the quality at least as good as Qp.
% We will use the matlab function quantile for it which takes two
% arguments, the vector with the results values and the probability of the
% quantile that you want to extract the value from.
% As a vector we take the final time that it has taken to converge if it is
% applicable. 
% After reading about quantiles a bit this seem useful but not just as a
% performance measure but to calculate probabilities of geting a value and
% extrapolate this to more trials and hence increasing the probability as
% defined before in the paper section. Nevertheless, this seems going out
% of track from the project so we will stick into calculate the probability
% of the quantile that has a maximum number of function evaluations as an
% indicative of probability of convergence in a different way (probability 
% of convergence in a desired number of function evaluations). 
% This can be used to analyse how much time does an algorithm needs to
% converge and how flexible this is (for some values of function evaluation
% the probability is going to be 0, none(NaN))

function [] = Quantiles (label, results, desiredT)

%%%%%%%%%%%%%%%%%%
% LAST CF VALUE  %
%%%%%%%%%%%%%%%%%%

vect = results(:,end).'; % Takes last column of results and converts it into a vector (just if necessary)

med = quantile(vect,0.5); % Takes the median of the dataset

disp(['The median is: ', num2str(med)]); 

p = linspace(0,1); % Linear space of probabilities from 0 to 1 taking 100 steps (default)
PDT =[]; % Probability(s) to the desired time 
for x=p   % For loop to calculate the probability associated to the quantile with our desired function evaluation value
    if round(quantile(vect,x)) == desiredT
        PDT = [PDT, x];
    end
end

PDT = sum(PDT)/length(PDT); % Average of the probabilities

disp(['The probability to converge with ', num2str(desiredT), ...
    ' function evaluations is: ', num2str(PDT)]);

mkdir(strcat('Quantiles',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Makes a new  folder that will contain the results
cd(strcat('Quantiles',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Changes directory to the new folder

save([label,'-Quantiles'],'PDT'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.

cd ..\; % Back to the original directory

end