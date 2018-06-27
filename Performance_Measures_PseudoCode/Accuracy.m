% Also known as relative error. It aims to determine where is located the 
% best solution found inside the interval defined by a lower bound 
% representing the worst known solution in the search space and an upper 
% bound representing the best known solution. Its advantadge is that
% minimizes posible biases due to normalization. The higher the obtained
% value the better performance of the algorithm. 
% For a time t the acuracy is the best solution in the interval minus the
% worst known solution divided by the best solution in the landscape minus
% the worst known solution. Again this will depend on the actual results
% that AMIGO gives, but the most probable is that the worst known solution
% is the starting point with the highest cost value and the best one is the
% value that we consider converged. As solution, it seems that we will have
% to use the cost value at each interval. 

function [] = Accuracy (label, results, worst, best) % Worst --> number of the column with the initial values of the cost function; Best --> value of the cost function at convergence(This value is sumed to be introduced by the user but it can aso be the index of the last column and there get the minimum value which will be the convergence value)

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

% Plot representing the accuracy over time (test)
x = [1:length(At)];
y = At;
figure
plot(x,y,'c');
title('Accuracy Plot')
xlabel('time')
ylabel('Accuracy')


mkdir(strcat('Accuracy',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Makes a new  folder that will contain the results
cd(strcat('Accuracy',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Changes directory to the new folder

save([label,'-Accuracy'],'At'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Acuracy over time is ', num2str(At)]); % This just prints the result on the screan to check it

cd ..\; % Back to the original directory


end