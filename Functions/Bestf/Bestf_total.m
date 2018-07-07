function [] = Bestf_total (directory)

cd ..\;

Bestf(directory);

cd ('Bestf');

SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Best_f')) 
        list2=[list2; SN(x,:)]; 
    end
end

Bf = [];
for x=1:length(list2(:,1))
    load(list2(x,:));
    Bf = [Bf,m];
end
BF = min(Bf);

save('BestTotal_f_Global', 'm');

disp(['The best cost function value (mnimum) is: ',num2str(BF)]);

end