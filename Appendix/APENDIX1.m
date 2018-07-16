%% APPENDIX 1: Simulations
%
% This appendix contains all the scripts used to run the different OID
% simulations for this work with comments on the function of the different
% parts of the scripts.

%% MPLac Model

% The script specifies the mathematical model structure as a charmodelC
% indicating name and number of states variables, model parameters,
% stimulus, observables, equations and parameter estimate values obtained
% by model fitting.

% Model introduction -- 'charmodelC'|'c_model'|'charmodelM'|
% 'matlabmodel'|'sbmlmodel'|'blackboxmodel'|'blackboxcost                        
model.input_model_type='charmodelC';                                             

% Number of states
model.n_st=3;                                                                     

% Number of model parameters 
model.n_par=8;                                                              

% Number of inputs, stimuli or control variables
model.n_stimulus=1;                                                            

% Names of the states
model.st_names=char('Cit_mrna','Cit_foldedP','Cit_fluo');                                                           

% Names of the parameters   
model.par_names=char('alpha1','Vm1','h1','Km1','d1','alpha2','d2','Kf');                      

% Names of the stimulus, inputs or controls 
model.stimulus_names=char('IPTG');                                                               

% Equations describing system dynamics. Time derivatives are regarded 'd'st_name''
model.eqns=...                                                              
               char('dCit_mrna=alpha1+Vm1*(IPTG^h1/(Km1^h1+IPTG^h1))-d1*Cit_mrna',...
                    'dCit_foldedP=alpha2*Cit_mrna-(d2+Kf)*Cit_foldedP',...
                    'dCit_fluo=Kf*Cit_foldedP-d2*Cit_fluo');           

% The parameter vector is set to the best estimates for MIP,r
model.par = [0.0164186333380725 0.291556643109224 1.71763487775568 5.14394334860864 0.229999999999978 6.63776658557266 0.00575139649497780 0.0216999999961899]; 

%% Steady State of the Model

% Function that calculates the steady state for the state variables (res)
% given a determined IPTG concentration and parameter values (theta) 
% introduced as inputs in the function

function [ res ] = InduciblePromoter_steady_state( theta, IPTG )

% Assignment of the different parameter values to their variables

a1 = theta(1);
Vm1 = theta(2);
h1 = theta(3);
Km1 = theta(4);
d1 = theta(5);
a2 = theta(6);
d2 = theta(7);
Kf = theta(8);

% Calculation of the steady states of the different variables setting the
% time derivative to 0.

Cit_mrna = (a1 + Vm1*(IPTG^h1/(Km1^h1+IPTG^h1)))/d1;

Cit_foldedP = (a2*Cit_mrna)/(Kf+d2);

Cit_fluo = (Kf*Cit_foldedP)/d2;

% Results included in a vector so it can be used outside the function.

res = [Cit_mrna Cit_foldedP Cit_fluo];

end

%% Initial Guesses Generator

% This function generates different initial guesses for the parameters in a
% logarithmic exploration from the maximum and minimum boundaries saving it
% in a matrix called ParFull. The number of vectors it is going to be the
% same as the number of experiments, in our case 30. 

function [] = MatrixParametersFunc (numExperiments)

% Selected boundaries for the parameters
theta_min = [3.88e-5,3.88e-2,0.5,2,7.7e-3,0.2433,5.98e-5,0.012];
theta_max = [0.4950,0.4950,4.9,10,0.23,6.8067,0.2449,0.0217];

% Create a matrix of initial guesses for the parameters, having as many
% rows as the number of PE iterations (numExperiments) 
% Each vector is passed as input to the computing function
M_norm = lhsdesign(numExperiments,length(theta_min));
M = zeros(size(M_norm));
for c=1:size(M_norm,2)
    for r=1:size(M_norm,1)
        M(r,c) = 10^(M_norm(r,c)*(log10(theta_max(1,c))-log10(theta_min(1,c)))+log10(theta_min(1,c))); % log exploration
    end
end 

ParFull = M; % in this case I am fitting all the values
save('MatrixParameters_InputComparison30.mat','ParFull');

end

%% Master Run Function

% This function runs the script run_in_silico_experiment_parfor_Optimised.m
% 20 times one after the other, each time with a different set of Global
% and Local optimisers and the different desired settings.

function [] = MasterRun2 ()

% Cell with the 20 different desired configurations of the inputs that 
% run_in_silico_experiment_parfor_Optimised.m requires in order. 

settings = {'1-eSS_dhc','eSS','dhc','','','',''; 
    '2-eSS_fminsearch','eSS','fminsearch','','','','';
    '3-eSS_nl2sol','eSS','nl2sol','','','','';
    '4-eSS_fmincon','eSS','fmincon','','','','';
    '5-DE_03056','de','',0.3,0.5,6,'';
    '6-DE_030856','de','',0.3,0.85,6,'';
    '7-DE_050856','de','',0.5,0.85,6,'';
    '8-DE_030853','de','',0.3,0.85,3,'';
    '9-DE_dhc_03056','de','dhc',0.3,0.5,6,'';
    '10-DE_nl2sol_030856','de','nl2sol',0.3,0.85,6,'';
    '11-DE_nl2sol_030852','de','nl2sol',0.3,0.85,2,'';
    '12-SRES_045','sres','','','','',0.45;
    '13-SRES_035','sres','','','','',0.35;
    '14-SRES_06','sres','','','','',0.6;
    '15-SRES_dhc_045','sres','dhc','','','',0.45;
    '16-SRES_nl2sol_045','sres','nl2sol','','','',0.45;
    '17-eSS_nolocal','eSS','','','','','';
    '18-DE_030859','de','',0.3,0.85,9,'';
    '19-DE_fminsearch_05056','de','fminsearch',0.5,0.5,6,'';
    '20-DE_fmincon_050856','de','fmincon',0.5,0.85,6,''};

% For loop that tries to run the script 
% run_in_silico_experiment_parfor_Optimised.m and if there is an error
% generates an error report and starts running the next simulation. 

for x=1:length(settings(:,1))
    try
        run_in_silico_experiment_parfor_Optimised(settings{x,1},settings{x,2},settings{x,3},settings{x,4},settings{x,5},settings{x,6},settings{x,7});
    catch err
        errorFile = ['SIMULATION-',x,'-ERROR.errorLog'];
        fid = fopen(errorFile,'a+');
        fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
        fclose(fid);
    end
end

end

%% Run Simulations Script

% Script that takes as an input the name or tab that will be included in
% the results name, the global solver selected, the local solver selected
% and the different settings for Differential Evolution (CR, F and
% strategy) and SRES (pf) to run the different simulations.

function [] = run_in_silico_experiment_parfor_Optimised(resultBase,globalSolv,localSolv,CR,F,strategy,pf)

% If statement to check that the Global solver name selected is correct
% since it is the most important part.
if strcmp(globalSolv, 'sres') || strcmp(globalSolv, 'de') || strcmp(globalSolv, 'eSS')
else
    fprintf('\nWrong Global Solver Selected. Please Try Again.\n\n');
    return
end

% Number of loops (1 since we are doing off-line OED) and the number of
% experiments performed for each simulation (30)
numLoops = 1; 
numExperiments = 30; 

% Start AMIGO2 Toolbox
AMIGO_Startup();

% Loads the Matrix with the different initial guesses for theta generated
% using the function MatrixParametersFunc.m
load('MatrixParameters_InputComparison30.mat');   

% Generates a new directory with the name introduced to resultBase and
% moves into it to save the results of the simulation in it.
mkdir(resultBase); 
cd(resultBase); 

% Parallel for loop (how many are done at the same time depends on the
% processor of the computer used) that runs a simulation 30 times (number
% of experiments)
parfor epcc_exps=1:numExperiments
        % Duration of the steps (200 min) and number of loops (1)
        stepd = 200;
        epccNumLoops = numLoops;
        
        % Try/catch statement to run the simulation from the script
        % fit_to_InduciblePromoter_Optimised_valuesOnly and if there is an
        % error and error file will be generated and the next trial will
        % start
        try
            global_theta_guess = ParFull(epcc_exps,:);
            % Generates the name of the file with the results
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

% Move back to the original folder
cd ..\; 

end

%% In silico experiment OID script 

% This script takes the same inputs that the function
% run_in_silico_experiment_parfor_Optimised plus the number of loops, step
% duration, number of experiments and name of the result file that are also
% present in the function. With this the script will realise OED to
% calculate the optimal input profile making use of the AMIGO Toolbox.

function [out]=fit_to_InduciblePromoter_Optimised_valuesOnly(epccOutputResultFileNameBase,epccNumLoops,stepd,epcc_exps,global_theta_guess,globalSolv,localSolv,CR,F,strategy,pf)

    % Settings for the output data file
    resultFileName = [strcat(epccOutputResultFileNameBase),'.dat'];
    rng shuffle;
    rngToGetSeed = rng;

    fid = fopen(resultFileName,'w');
    fprintf(fid,'HEADER DATE %s\n', datestr(datetime()));
    fprintf(fid,'HEADER RANDSEED %d\n',rngToGetSeed.Seed);
    fclose(fid);

    startTime = datenum(now);

    % Clear previous variables so it does not interfere    
    clear model;
    clear exps;
    clear best_global_theta_log;
    clear pe_results;
    clear ode_results;

    % Where results will be stored    
    results_folder = strcat('InduciblePromoter',datestr(now,'yyyy-mm-dd-HHMMSS'));
    short_name     = strcat('IndProm',int2str(epcc_exps));

    % Read the model into the model variable
    InduciblePromoter_load_model;

    % We start with no experiments
    exps.n_exp=0;

    % Defining boundaries for the parameters 
    global_theta_max = [0.4950,0.4950,4.9,10,0.23,6.8067,0.2449,0.0217];       
    global_theta_min = [3.88e-5,3.88e-2,0.5,2,7.7e-3,0.2433,5.98e-5,0.012];    
    global_theta_guess = global_theta_guess';

    % Specify the parameters to be calibrated.
    param_including_vector = [true,true,true,true,true,true,true,true];

    % Compile the model
    clear inputs;
    inputs.model = model;
    inputs.pathd.results_folder = results_folder;
    inputs.pathd.short_name     = short_name;
    inputs.pathd.runident       = 'initial_setup';
    AMIGO_Prep(inputs);

    % Duration of the experiment (3000s)
    totalDuration = 50*60; 
    % Number of loops (1)
    numLoops = epccNumLoops;           
    % Duration of each loop
    duration = totalDuration/numLoops;
    % Duration of each steep (200 min)
    stepDuration = stepd;                

for i=1:numLoops
    
    %Compute steady state 
    y0 = InduciblePromoter_steady_state(global_theta_guess,0);
        
    if exps.n_exp == 0
        oid_y0 = [y0]; 
        best_global_theta = global_theta_guess;
    else
        % Simulate the experiment without noise to find end state
        clear inputs;
        inputs.model = model;
        inputs.model.par = best_global_theta;
        inputs.exps = exps;
        
        inputs.plotd.plotlevel='noplot';
        inputs.pathd.results_folder = results_folder;
        inputs.pathd.short_name     = short_name;
        inputs.pathd.runident       = strcat('sim-',int2str(i));
        
        sim = AMIGO_SData(inputs);
        
        oid_y0 = [sim.sim.states{1}(end,:)];
        
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimal experiment design %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear inputs;
    inputs.model = model;
    inputs.exps  = exps;
    format long g
        
    % Number of experiments to be designed (1)
    inputs.exps.n_exp = inputs.exps.n_exp + 1; 
    % Index of experiment
    iexp = inputs.exps.n_exp;    
    % Type of experiment. od = optimaly designed
    inputs.exps.exp_type{iexp}='od'; 
    % Number  of observables
    inputs.exps.n_obs{iexp}=1;    
    % Name of observables
    inputs.exps.obs_names{iexp}=char('Fluorescence');   
    % Observables function
    inputs.exps.obs{iexp}=char('Fluorescence = Cit_fluo');          
        
    % Initial conditions
    inputs.exps.exp_y0{iexp}=oid_y0; 
    % Duration of experiment (3000 min)
    inputs.exps.t_f{iexp}=duration;
    % Number of sampling times (1 every 5 min)
    inputs.exps.n_s{iexp}=duration/5+1;                             
    
    % OED of the input
    inputs.exps.u_type{iexp}='od';
    % Stimuli is steps with fixed duration
    inputs.exps.u_interp{iexp}='stepf';   
    % Number of steps (15)
    inputs.exps.n_steps{iexp}=round(duration/stepDuration);  
    % Switching times
    inputs.exps.t_con{iexp}=0:stepDuration:(duration);    
    % Minimum value for the input
    inputs.exps.u_min{iexp}=0*ones(1,inputs.exps.n_steps{iexp});    
    % Maximum value for the input
    inputs.exps.u_max{iexp}=1000*ones(1,inputs.exps.n_steps{iexp}); 
    
    % Loading information about the parameters
    inputs.PEsol.id_global_theta=model.par_names(param_including_vector,:);
    inputs.PEsol.global_theta_guess=transpose(global_theta_guess(param_including_vector));
    inputs.PEsol.global_theta_max=global_theta_max(param_including_vector);  
    inputs.PEsol.global_theta_min=global_theta_min(param_including_vector);  
       
    % Deffinition of the exerimental noise (homoscedastic)
    inputs.exps.noise_type='homo_var';           
    inputs.exps.std_dev{iexp}=[0.05];
    inputs.OEDsol.OEDcost_type='Dopt';
    
    % SIMULATION
    % IVP solver -- 'cvodes'(default, C)|'ode15s' 
    % (default, MATLAB, sbml)|'ode113'|'ode45'
    inputs.ivpsol.ivpsolver='cvodes';  
    % Sensitivities solver: 'cvodes'(default, C)| 'sensmat'(matlab)|
    % 'fdsens2'|'fdsens5'
    inputs.ivpsol.senssolver='cvodes';                    
    % IVP solver integration tolerances
    inputs.ivpsol.rtol=1.0D-8;                            
    inputs.ivpsol.atol=1.0D-8;
    
    % OPTIMIZATION
    
    % If statement to select the settings for the Global Optimiser selected
    % If eSS is selected then the next bit of code is considered
    if strcmp(globalSolv, 'eSS')
        % Deffinition of the global solver (eSS)
        inputs.nlpsol.nlpsolver=globalSolv;
        % Maximum number of function evaluations
        inputs.nlpsol.eSS.maxeval = 3e5;
        % Maximum number of CPU time
        inputs.nlpsol.eSS.maxtime = 30e3;
        % Definition of the local solver for the hybridisation
        inputs.nlpsol.eSS.local.solver = localSolv; 
        inputs.nlpsol.eSS.local.finish = localSolv;    
        % If nl2sol is selected as a local solver the next bit of code is
        % considered
        if strcmp(localSolv, 'nl2sol')
            % Maximum number of iterations and function evaluations for
            % nl2sol
            inputs.nlpsol.eSS.local.nl2sol.maxiter  =     500;     
            inputs.nlpsol.eSS.local.nl2sol.maxfeval =     500;     
        end
        inputs.nlpsol.eSS.log_var=1:inputs.exps.n_steps{iexp};
    % If de is selected then the next bit of code is considered
    elseif strcmp(globalSolv, 'de') 
        % If a loca solver is selected then hybridisation is performed
        if strcmp(localSolv, 'dhc') || strcmp(localSolv, 'n12sol')|| strcmp(localSolv, 'fminsearch')|| strcmp(localSolv, 'fmincon')
            hyb = ['hyb_de_',localSolv];
            inputs.nlpsol.nlpsolver=hyb;
        else    
            inputs.nlpsol.nlpsolver=globalSolv;
        end
        % Population number size (usually greater than 10*number of 
        % decision variables)
        inputs.nlpsol.DE.NP = max([200, 10*(2*inputs.exps.n_steps{iexp}-1)]);
        % Maximum number of iterations (number of funtion evaluations is
        % population size times the number of iterations
        inputs.nlpsol.DE.itermax = round((500*1e3)/inputs.nlpsol.DE.NP);   
        % Maximum variance of population
        inputs.nlpsol.DE.cvarmax = 1e-5;
        % DE_Stepsize
        inputs.nlpsol.DE.F = F;                               
        % Crossover constant
        inputs.nlpsol.DE.CR = CR;      
        % Strategy selected.
        % 1 --> DE/best/1/exp
        % 2 --> DE/rand/1/exp
        % 3 --> DE/rand-to-best/1/exp
        % 4 --> DE/best/2/exp
        % 5 --> DE/rand/2/exp
        % 6 --> DE/best/1/bin
        % 7 --> DE/rand/1/bin
        % 8 --> DE/rand-to-best/1/bin
        % 9 --> DE/best/2/bin
        % else DE/rand/2/bin
        inputs.nlpsol.DE.strategy = strategy;
        inputs.nlpsol.DE.refresh=2; 
    % If sres is selected then the next bit of code is selected
    elseif strcmp(globalSolv, 'sres') 
        % If a loca solver is selected then hybridisation is performed
        if strcmp(localSolv, 'dhc') || strcmp(localSolv, 'n12sol')|| strcmp(localSolv, 'fminsearch')|| strcmp(localSolv, 'fmincon')
            hyb = ['hyb_sres_',localSolv];
            inputs.nlpsol.nlpsolver=hyb;
        else    
            inputs.nlpsol.nlpsolver=globalSolv;
        end
        % Population size
        inputs.nlpsol.SRES.NP = max([200, 10*(2*inputs.exps.n_steps{iexp}-1)]);
        % Maximum number of iterations (number of funtion evaluations is
        % population size times the number of iterations
        inputs.nlpsol.SRES.itermax = round((500*1e3)/inputs.nlpsol.SRES.NP);     
        % Maximum variance of population
        inputs.nlpsol.SRES.cvarmax = 1e-5;        
        % Expected rate of convergence
        inputs.nlpsol.SRES.varphi = 1;                          
        % Pressure of fitness
        inputs.nlpsol.SRES.pf = pf;
    end
    
    % Specifications for the results
    inputs.plotd.plotlevel='noplot';
    
    inputs.pathd.results_folder = results_folder;
    inputs.pathd.short_name     = short_name;
    inputs.pathd.runident       = strcat('oed-',int2str(i));
        
    oed_start = now;
    
    results = AMIGO_OED(inputs);
    oed_results{i} = results;
    oed_end = now;
    
    results.plotd.plotlevel = 'noplot';
    
end   

%Saving of the results
true_param_values = model.par(param_including_vector);

save([strcat(epccOutputResultFileNameBase),'.mat'], 'oed_results','exps','inputs','true_param_values','best_global_theta','best_global_theta_log');

out= 1;
end