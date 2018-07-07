function [] = Bestf (directory)

cd ..\;
cd (directory);

SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end


f_val = [];
for x=1:length(list2(:,1))
    load(list2(x,:));
    f = oed_results{1}.nlpsol.fbest;
    f_val = [f_val, f];
end

m = min(f_val);

cd ..\;
cd ('Functions');
cd ('Bestf');
save(['Best_f_',directory], 'm');
cd ..\

end