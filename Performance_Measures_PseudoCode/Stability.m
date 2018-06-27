% Reflection of how much an algorithm is able to recover after a change. 
% It is accuracy dependent since an algorithm is said to be stable if it 
% maintains its accuracy from a time step to the next one. The most stable 
% algorithm is not necessary the best performing algorithm in terms of 
% producing quality solutions. 
% At a generation g is the difference between the current accuracy and the
% accuracy at the previous generation(g-1) if this difference is positive
% or 0 otherwise. The closer to zero, the higher is the stability.

function [] = Stability (label, results, worst, best)

%%%%%%%%%%%%%%
%  ACCURACY  %
%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%
%  STABILITY  %
%%%%%%%%%%%%%%%

Sg = []; % Empty vector that will contain all the values for the stability over generations

for x=1:length(At)
    if x-1 ~= 0 % To avoid error with the first value which will not have a previous one
        s = At(x)-At(x-1); % Calculation of the stability values at each generation
        if s < 0 % This bit of code substitutes all negative values for zero as the definition of stability indicates
            s = 0;
        end
        Sg = [Sg, s]; % Addition of all the different values into the vector Sg
    end
end

% Plot representing the accuracy over time (test)
x = [1:length(Sg)];
y = Sg;
figure
plot(x,y,'c');
title('Stability Plot')
xlabel('generation')
ylabel('Stability')


mkdir(strcat('Stability',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Makes a new  folder that will contain the results
cd(strcat('Stability',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Changes directory to the new folder

save([label,'-Stability'],'Sg'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Stability over generations is ', num2str(Sg)]); % This just prints the result on the screan to check it

cd ..\; % Back to the original directory


end
