% radarplot([1 ; 2; 0 ; 3; 2 ],{ 'a' ; 'b' ; 'c' ; 'd' ; 'e' } ,'r')


function [] = RadarPlots (label)
cd ('R-Measure');
SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),label)) 
        list2=[list2; SN(x,:)]; 
    end
end
for x=1:length(list2(:,1))
    load(list2(x,:));
end
cd ..\;
clear list2;

cd ('Q_Measure');
SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),label)) 
        list2=[list2; SN(x,:)]; 
    end
end
for x=1:length(list2(:,1))
    load(list2(x,:));
end
clear list2;
cd ..\;

cd ('MeanFunctionEvaluations');
SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),label)) 
        list2=[list2; SN(x,:)]; 
    end
end
for x=1:length(list2(:,1))
    load(list2(x,:));
end
clear list2;
cd ..\;

meas1 = [msd;sd;Qm];
tags1 = {'Modified Standard Deviation';'Standard Deviation';'Q Measure'};
meas2 = [mvr;vr;Qm];
tags2 = {'Modified Variance';'Variance';'Q Measure'};

radarplot(meas1,tags1 ,'r');
title(label);

radarplot(meas2,tags2 ,'r');
title(label)

end