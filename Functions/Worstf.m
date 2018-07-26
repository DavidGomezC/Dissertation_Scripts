function [] = Worstf (directory)

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
f_val_De = [];
for x=1:length(list2(:,1))
    if (contains(list2(x,:),'eSS'))
        load(list2(x,:));
        f = oed_results{1}.nlpsol.f(1);
        f_val = [f_val, f];
    elseif (contains(list2(x,:),'DE'))
        load(list2(x,:));
        f = oed_results{1}.nlpsol.conv_curve(:,2);
        f_val_De = [f_val_De, f];
    end
end

if ~isempty(f_val)
    w = max(f_val);
elseif ~isempty(f_val_De)
    w = max(f_val_De(1,:));
end

cd ..\;
cd ('Functions');
cd ('Worstf');
save(['Worst_f_',directory], 'w');
cd ..\

end
