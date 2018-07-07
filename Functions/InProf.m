function [] = InProf (directory)

cd ..\;
cd (directory);

best = [];
SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

oip = {};
roip= {};
for x=1:length(list2(:,1))
    load(list2(x,:));
    f = oed_results{1}.nlpsol.fbest;
    best = [best, f];
    us = oed_results{1}.oed.u{1};
    roip = [roip,us];
    u=[];
    for x=1:length(us)
        r = repelem(us(x),200);
        u = [u, r];
    end
    oip{end+1} = u;
end


cd ..\;
cd ('Functions');

cd ('InputProfiles');
for x=1:length(oip)
    fig = plot(1:length(oip{x}), oip{x});
    title(['Input Profile ', num2str(x), ' Cost Function: ', num2str(best(x))]);
    xlabel('Time (s)');
    ylabel('[IPTG]');
    saveas(fig,[num2str(x),'_feval.png']);
end
cd ..\
end
