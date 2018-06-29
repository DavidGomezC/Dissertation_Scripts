% Quantifies the distance between the performance curves of each pair of
% algorithms instead of measuring their closeness to the optimum. The area
% between curves for two algorithms is the integral of the difference of
% two functions that are a measure of population qualty. This measure qill
% be the accuracy as the two definitions given in fitness degradation (the 
% second one can be used since it is similar to the average best of 
% generation).
% Since we don't have the function that descrives them but the actual
% values at each point we can use the function trapz from matlab which will
% aproximate the area under the curves as a numeric value and then we can
% substract this two to obtain the area between curves.
% A positive value means that the area under the curve of the first
% algorithm is higher, a negative value means that the second is higher. 
% This can be used if necessary to compare two algorithms when the usual
% measures cannot distinguish between them. 
% A plot with the area between curves colored will be also done for a
% visual representation of the concept and results


function [] = AreaBetweenCurves (label, resultsA1, resultsA2) % resultsA1 is the vector with the accuracy results for the frst algorithm and resultsA2 for the second algorithm

x1 = 1:length(resultsA1); % Define the x axis of the function and plot
x2 = 1:length(resultsA2);

pA1 = trapz(x1,resultsA1); % Calculation of the area under each curve
pA2 = trapz(x2,resultsA2);


ABC = pA1 - pA2; % Calculation of the area between curves

plot(x1, resultsA1, 'LineWidth',3); % Plot to visualize the results
hold on;
plot(x2, resultsA2, 'LineWidth',3);
patch([x1 fliplr(x1)], [resultsA1 fliplr(resultsA2)], 'g');
title(['Area between curves: ', num2str(ABC)]);
xlabel('generation');
ylabel('Accuracy');
hold off;


mkdir(strcat('AreaBetweenCurves',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Makes a new  folder that will contain the results
cd(strcat('AreaBetweenCurves',datestr(now,'yyyy-mm-dd-HHMMSS'))); % Changes directory to the new folder

save([label,'-AreaBetweenCurves'],'ABC'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['The area between curves is ', num2str(ABC)]); % This just prints the result on the screan to check it


cd ..\; % Back to the original directory


end