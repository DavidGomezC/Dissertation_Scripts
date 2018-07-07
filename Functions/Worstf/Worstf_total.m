function [] = Worstf_total (directory)

cd ..\;

Worstf(directory);

cd ('Worstf');

SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Worst_f')) 
        list2=[list2; SN(x,:)]; 
    end
end

Wf = [];
for x=1:length(list2(:,1))
    load(list2(x,:));
    Wf = [Wf,w];
end
WF = min(Wf);

save('WorstTotal_f_Global', 'w');

disp(['The Worst cost function value (mnimum) is: ',num2str(WF)]);

end