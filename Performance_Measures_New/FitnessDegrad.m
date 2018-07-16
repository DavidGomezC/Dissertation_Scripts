function [] = FitnessDegrad (directory)

%%%%%%%%%%%%%%%
%  ACCURACY   %
%%%%%%%%%%%%%%%

Accuracy(directory);
cd ('Accuracy\');
load([directory,'-Accuracy']);
cd ..\;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR REGRESSION ACCURACY   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = [1:length(At)];
y = At;

coef = polyfit(x,y,1);
rl1 = coef(1)*x + coef(2);
beta = coef(1);
figure;
hold on
fig1 = plot(x,y);
fig1 = plot(x,rl1);
title(['Linear Regression 1. Beta = ', num2str(beta)]);
xlabel('Periods');
ylabel('Accuracy');
hold off
save([directory,'-Beta'],'beta'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Beta degradation is ', num2str(beta)]); % This just prints the result on the screan to check it
saveas(fig1,[directory,'_BDegr1.png']);

cd ..\; % Back to the original directory

%%%%% Save Results CSV %%%%%
A1 = [1, y].';
R1 = [2, rl1].';
FD = [A1(:), R1(:)];
csvwrite([directory,'FitnessDegradation.csv'],FD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
