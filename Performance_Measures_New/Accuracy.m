function [] = Accuracy (directory) % Worst --> number of the column with the initial values of the cost function; Best --> value of the cost function at convergence(This value is sumed to be introduced by the user but it can aso be the index of the last column and there get the minimum value which will be the convergence value)
% ACCESS BEST SOLUTION
cd('..\Functions\Bestf');
load('BestTotal_f_Global');
bsl = m;

% ACCESS WORST SOLUTION
cd ('..\Worstf');
load('WorstTotal_f_Global');
wks = w;

% GET THE COST FUNCTIONS OF INTEREST
cd ('..\..\');
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

BS = min(F);

At = (BS-wks)/(bsl-wks);


% Plot representing the accuracy over time (test)
x = 1:length(At);
y = At;
figure
fig = plot(x,y);
title('Accuracy Plot')
xlabel('Function Evaluations')
ylabel('Accuracy')

cd(strcat('Accuracy')); % Changes directory to the new folder

save([directory,'-Accuracy'],'At'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.
disp(['Acuracy al the last time is time is ', num2str(At(end))]); % This just prints the result on the screan to check it
saveas(fig,[directory,'_feval.png']);

cd ..\; % Back to the original directory

end