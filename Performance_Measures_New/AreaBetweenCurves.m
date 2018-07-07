function [] = AreaBetweenCurves (resultsA1, resultsA2) % resultsA1 is the vector with the accuracy results for the frst algorithm and resultsA2 for the second algorithm

cd('Accuracy');

SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:), resultsA1)) || (contains(SN(x,:),resultsA2))
        list2=[list2; SN(x,:)]; 
    end
end

curvs = {};
for x = 1:length(list2(:,1))
    load(list2(x,:));
    curvs{x} = At;
end

cd ..\;
x1 = 1:length(curvs{1}); % Define the x axis of the function and plot
x2 = 1:length(curvs{2});
y1 = curvs{1};
y2 = curvs{2};

pA1 = trapz(x1,y1); % Calculation of the area under each curve
pA2 = trapz(x2,y2);
 
ABC = pA1 - pA2; % Calculation of the area between curves

% hold on;
% plot(x1, y1, 'LineWidth',3); % Plot to visualize the results
% 
% plot(x2, y2, 'LineWidth',3);
% patch([x1 fliplr(x2)], [y1 fliplr(y2)], 'g');
% title(['Area between curves: ', num2str(ABC)]);
% xlabel('generation');
% ylabel('Accuracy');
% hold off;

cd(strcat('AreaBetweenCurves')); % Changes directory to the new folder

save([resultsA1,'-', resultsA2,'-AreaBetweenCurves'],'ABC'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['The area between curves is ', num2str(ABC)]); % This just prints the result on the screan to check it


cd ..\; % Back to the original directory


end