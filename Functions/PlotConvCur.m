function [] = PlotConvCur (directory)

cd ..\;
cd (directory);

SN = ls; 
list2 =[]; 
for x=1:length(SN(:,1)) 
    if (contains(SN(x,:),'Optsteps')) 
        list2=[list2; SN(x,:)]; 
    end
end

figure1 = {};
figure2 = {};
for x=1:length(list2(:,1))
    load(list2(x,:));
    xeval = oed_results{1}.nlpsol.neval;
    xtime = oed_results{1}.nlpsol.time;
    y = oed_results{1}.nlpsol.f;
    figure1{end+1} = xeval;
    figure1{end+1} = y;
    figure2{end+1} = xtime;
    figure2{end+1} = y;
end

figure
hold on
for x=1:2:length(figure1)
    fig1 = plot(figure1{x},figure1{x+1});
end
title('Convergence curve');
xlabel('Function evaluations');
ylabel('f');
hold off


figure
hold on
for x=1:2:length(figure2)
    fig2 = plot(figure2{x},figure2{x+1});
end
title('Convergence curve');
xlabel('CPU time');
ylabel('f');
hold off


figure
hold on
for x=1:2:length(figure1)
    fig3 = plot(figure1{x},figure1{x+1});
end
set(gca, 'YScale', 'log');
title('Convergence curve');
xlabel('Function evaluations');
ylabel('f');
hold off

cd ..\;
cd ('Functions');

cd ('Convergence_Curves_Fig');
saveas(fig1,[directory,'_feval.png']);
saveas(fig2,[directory,'_time.png']);
saveas(fig3,[directory,'_feval_log.png']);
cd ..\

end