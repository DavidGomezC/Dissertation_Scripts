function [] = FitnessDegrad (directory)

%%%%%%%%%%%%%%%%
%  ACCURACY 1  %
%%%%%%%%%%%%%%%%

Accuracy(directory);
cd ('Accuracy\');
load([directory,'-Accuracy']);
cd ..\;

%%%%%%%%%%%%%%%%
%  ACCURACY 2  %
%%%%%%%%%%%%%%%%

% GET THE COST FUNCTIONS OF INTEREST
cd ('..\');
cd (directory);

SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

cf = {};
fe={};
for x=1:length(list2(:,1))
    load(list2(x,:));
    f = oed_results{1}.nlpsol.f;
    s = oed_results{1}.nlpsol.neval;
    cf{x} = f;
    fe{x} = s;
end

cd ..\; 
cd ('Performance_Measures_New');

x = 1:300800;
cft = {};

parfor q=1:length(cf)
    ind = [];
    for w=1:length(x)
        if ismember(x(w), fe{q})
            ind = find(fe{q}==x(w));
            cft{q}(w) = cf{q}(ind);
        else
            cft{q}(w) = nan;
        end
    end    
end

F = [];
for q=1:length(cft)
    F = [F; fillmissing(cft{q},'previous')]; % It can be previous or linear
end
F(isnan(F))=0; 
F(2,:)=[];
Ak = (sum(F))/(length(F(:,1))); % Empty vector that will have the cost function values but just the best ones for each generation
 
cd(strcat('Fitness_Degradation')); % Changes directory to the new folder

save([directory,'-Acuracy'],'At'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
save([directory,'-Acuracy2'],'Ak'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR REGRESSION ACCURACY 1  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
save([directory,'-Beta1'],'beta'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Beta degradation 1 is ', num2str(beta)]); % This just prints the result on the screan to check it
saveas(fig1,[directory,'_BDegr1.png']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LINEAR REGRESSION ACCURACY 2  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = [1:length(Ak)];
t = Ak;

coef2 = polyfit(x,t,1);
rl2 = coef2(1)*x + coef2(2);
beta2 = coef2(1);
figure;
hold on
fig2 = plot(x,t);
fig2 = plot(x,rl2);
title(['Linear Regression 2. Beta = ', num2str(beta2)]);
xlabel('Periods');
ylabel('Accuracy');
hold off

save([directory,'-Beta2'],'beta2'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Beta degradation 2 is ', num2str(beta2)]); % This just prints the result on the screan to check it
saveas(fig2,[directory,'_BDegr2.png']);

cd ..\; % Back to the original directory

%%%%% Save Results CSV %%%%%
A1 = [1, y].';
R1 = [2, rl1].';
A2 = [3, t].';
R2 = [4, rl2].';
FD = [A1(:), R1(:), A2(:), R2(:)];
csvwrite([directory,'FitnessDegradation.csv'],FD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
