function [] = InputProf ()

best = -5e+42;

SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end


for x=1:length(list2(:,1))
    load(list2(x,:));
    if oed_results{1}.nlpsol.fbest <= best
        IP = oed_results{1}.nlpsol.vbest;
    end
end
I = [0, IP];
disp(I);
u = [];
for x=1:length(IP)
    r = repelem(IP(x),200);
    u = [u, r];
end
csvwrite('InputSing.csv',I.');
% tag = 'InputP';
%IP = [tag, IP];
disp(length(u.'));






end