function [] = run_in_silico_experiment_parfor_Optimised(resultBase,globalSolv,localSolv,CR,F,strategy,pf)

if strcmp(globalSolv, 'sres') || strcmp(globalSolv, 'de') || strcmp(globalSolv, 'eSS')
else
    fprintf('\nWrong Global Solver Selected. Please Try Again.\n\n');
    return
end

numLoops = 1; % Number of loops fixed inside the function instead of an input because it will allways be the same
numExperiments = 30; % Number of experiments fixed inside the function instead of an input because it will allways be the same

AMIGO_Startup();

load('MatrixParameters_InputComparison30.mat');   

mkdir(resultBase); % New folder where to keep all the .dat files instead of in the same folder and with the tag name used for the simulation for track easines
cd(resultBase); % Change directory to the neww folder

parfor epcc_exps=1:numExperiments
        stepd = 200;
        epccNumLoops = numLoops;
        try
            global_theta_guess = ParFull(epcc_exps,:);
            epccOutputResultFileNameBase = [resultBase,'-','Optsteps_',globalSolv,'_',localSolv,'-',num2str(numLoops),'_loops-',num2str(epcc_exps)];
            [out]=fit_to_InduciblePromoter_Optimised_valuesOnly(epccOutputResultFileNameBase,epccNumLoops,stepd,epcc_exps,global_theta_guess,globalSolv,localSolv,CR,F,strategy,pf);
            
        catch err
            %open file
            errorFile = [resultBase,'-','Optsteps_',globalSolv,'-','_',localSolv,num2str(numLoops),'_loops-',num2str(epcc_exps),'.errorLog'];
            fid = fopen(errorFile,'a+');
            fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
            % close file
            fclose(fid);
        end

end

% Lines of code necessary to set up everything to be able to send emails
% from the scripts and aboid errors
setpref('Internet','E_mail','david.andorra93@gmail.com');
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username','david.andorra93@gmail.com');
setpref('Internet','SMTP_Password','pasword');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
props.setProperty( 'mail.smtp.starttls.enable', 'true' );

% This bit of codde is set to send an email. If everything goes well, it
% will send an email saying that the simulation has finished, but if the
% simulation stops early due to errors and generates a errorLog file the
% mail will indicate it so we know that it has finished due to errors, no a
% succesful simulation
sn = ls;
for x=1:length(sn(:,1))
    if contains(sn(x,:),'.errorLog')
        sendmail('s1778490@sms.ed.ac.uk','SOMETHING FAILED!', ...
            'Some simulation has failed generating an errorLog file');
        break
    else
        sendmail('s1778490@sms.ed.ac.uk','SIMULATION FINISHED!', ...
            ['The simulation ', resultBase ,' has finished!']);
        break
    end
end

%cd ..\; % Move back to the original folder

end
