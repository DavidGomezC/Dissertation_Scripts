%% Run Quantiles

% This function runs the functions Quantiles.m and Quantiles2.m using 8
% different values at each time. As inputs it takes a the name of the 
% directory where the AMIGO2 results desired to analyse are.

function [] = RunQuantiles (directory)

% Heather
fprintf('\nQuantiles1: \n\n');
% Linear space of different cost function exponents to calculate its
% probability using the function Quantiles.m stored in the variable q1.
q1 = linspace(18,42,8);
% For loop calling the function Quantiles.m for each value stored in q1.
for x = q1
    Quantiles(directory,round(x));
end

% Heather
fprintf('\nQuantiles2: \n\n');
% Liner space of different function evaluation values to calculate its
% probability of convergence using the function Quantiles2.m stored in the
% variable q2.
q2 = linspace(800,200100,8);
% For loop calling the function Quantiles2.m for each value stored in q2.
for x = q2
    Quantiles2(directory,round(x));
end
end