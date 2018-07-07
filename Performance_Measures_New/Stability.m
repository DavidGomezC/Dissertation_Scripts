function [] = Stability (directory)

%%%%%%%%%%%%%%
%  ACCURACY  %
%%%%%%%%%%%%%%

Accuracy(directory);
cd ('Accuracy\');
load([directory,'-Accuracy']);
cd ..\;

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
fig = plot(x,y,'c');
title('Stability Plot')
xlabel('Function Evaluation')
ylabel('Stability')


cd(strcat('Stability')); % Changes directory to the new folder

save([directory,'-Stability'],'Sg'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
saveas(fig,[directory,'_feval.png']);

cd ..\; % Back to the original directory


end
