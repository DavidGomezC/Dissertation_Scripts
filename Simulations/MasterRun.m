function [] = MasterRun ()

run_in_silico_experiment_parfor_Optimised('eSS_dhc','eSS','dhc','','','','');

run_in_silico_experiment_parfor_Optimised('eSS_fminsearch','eSS','fminsearch','','','','');

run_in_silico_experiment_parfor_Optimised('eSS_nl2sol','eSS','nl2sol','','','','');

run_in_silico_experiment_parfor_Optimised('eSS_fmincon','eSS','fmincon','','','','');

run_in_silico_experiment_parfor_Optimised('DE_03056','de','',0.3,0.5,6,'');

run_in_silico_experiment_parfor_Optimised('DE_030856','de','',0.3,0.85,6,'');

run_in_silico_experiment_parfor_Optimised('DE_050856','de','',0.5,0.85,6,'');

run_in_silico_experiment_parfor_Optimised('DE_030853','de','',0.3,0.85,3,'');

run_in_silico_experiment_parfor_Optimised('DE_dhc_03056','de','dhc',0.3,0.5,6,'');

run_in_silico_experiment_parfor_Optimised('DE_nl2sol_030856','de','nl2sol',0.3,0.85,6,'');

run_in_silico_experiment_parfor_Optimised('DE_nl2sol_030852','de','nl2sol',0.3,0.85,2,'');

run_in_silico_experiment_parfor_Optimised('SRES_045','sres','','','','',0.45);

run_in_silico_experiment_parfor_Optimised('SRES_035','sres','','','','',0.35);

run_in_silico_experiment_parfor_Optimised('SRES_06','sres','','','','',0.6);

run_in_silico_experiment_parfor_Optimised('SRES_dhc_045','sres','dhc','','','',0.45);

run_in_silico_experiment_parfor_Optimised('SRES_nl2sol_045','sres','nl2sol','','','',0.45);

run_in_silico_experiment_parfor_Optimised('eSS_nolocal','eSS','','','','','');

run_in_silico_experiment_parfor_Optimised('DE_030859','de','',0.3,0.85,9,'');

run_in_silico_experiment_parfor_Optimised('DE_fminsearch_05056','de','fminsearch',0.5,0.5,6,'');

run_in_silico_experiment_parfor_Optimised('DE_fmincon_050856','de','fmincon',0.5,0.85,6,'');

end