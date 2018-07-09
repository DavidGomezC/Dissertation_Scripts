function [] = ConvCurvToCSV (label)

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

for x=1:length(cf)
    cc = [];
    cc = [fe{x}; cf{x}];
    disp(cc);
    csvwrite([label, num2str(x),'.csv'],cc);

end


end