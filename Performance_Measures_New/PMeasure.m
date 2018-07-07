% Hypothesis --> We can modify it as a mesure of convergence rate from the 
% point of view of each stocstic iteration. Population will refear to the number of trials then. 
% So, instead of the average
% vector of individuals we can use the average cost value per function
% evaluation, so instead of per population members we consider per trials.
% The center of the population will be average cost function value for the
% total amount of trials. Hence P-measure will be the euclidean distance
% between the center of a cost value population at each function evaluation
% and the individual fardest from it. 

function [] = PMeasure (directory, label)

ind={}; % Empty cell of individuals

cd ..\;
cd (directory);

SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

cf = {};
fe = {};
for x=1:length(list2(:,1))
    load(list2(x,:));
    s = oed_results{1}.nlpsol.neval;
    c = oed_results{1}.nlpsol.f;
    cf{x} = c;
    fe{x} = s;
end

cd ..\;
cd ('Performance_Measures_New');

x = 1:300800;
y = {};

parfor q=1:length(cf)
    ind = [];
    for w=1:length(x)
        if ismember(x(w), fe{q})
            ind = find(fe{q}==x(w));
            y{q}(w) = cf{q}(ind);
        else
            y{q}(w) = nan;
        end
    end    
end

F = [];
for q=1:length(y)
    F = [F; fillmissing(y{q},'previous')]; % It can be previous or linear
end
F(isnan(F))=0; 

OP = sum(F)/length(F(:,1));

% Op=[]; % Empty vector with the center of the population
% parfor q=1:length(x) % For loop that goes through the length of the vector of individuals
%     o=sum(F(:,q))/length(F(:,q));
%     Op = [Op, o];
% end

dis = abs(F - OP); 

Pm = max(dis); 


cd(strcat('P-Measure')); % Changes directory to the new folder

save([label,'-PMeasure'],'Pm'); % This saves the value of the probability calculated in the current folder as a .mat file and with the label at the beginning to be able to identify it.

cd ..\; % Back to the original directory
plot(Pm);
disp('Everything correct!');

end
