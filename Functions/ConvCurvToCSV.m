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
    if (contains(list2(x,:),'eSS'))
        load(list2(x,:));
        s = oed_results{1}.nlpsol.neval;
        c = oed_results{1}.nlpsol.f;
        cf{x} = c; fe{x} = s;
    elseif (contains(list2(x,:),'DE'))
        load('FunctionEvaluationsVector');
        load(list2(x,:));       
        s = FE;
        c = oed_results{1}.nlpsol.conv_curve(:,2);
        cf{x} = c.'; fe{x} = s;
    end
end

for x=1:length(cf)
    cc = [];
    cc = [fe{x}; cf{x}];
    disp(cc);
    csvwrite([label, num2str(x),'.csv'],cc);

end


end
