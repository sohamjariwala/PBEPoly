% Script to optimize the values for parameters in the polydisperse
% population balance consitutive model.
clear; clc
Nfits = 15;

% Initial vector of parameters
numParameters = 4;

% Initial guess
parVec =[0.016293920207754  0.6 -6.863277109520129   0.857693893803876];

lb = min([0.75*parVec; 1.25*parVec]);
ub = max([0.75*parVec; 1.25*parVec]);

% Fitting using Genetic algorithm
W = [];
alfa = W; b_0 = W; d_f = W;
porosity = W; m_p = W;
error = W; N = W; G_0 = W; sigma_y0 = W; mu_s = W;

% i =1;
for i = 1:Nfits
fprintf("Starting Run %d\n", i);

% Initializing the fluid object and functions with default parameter
fluid(i) = PBEPoly;
% Changing constants and parameters to ones obtained from monodisperse 
% solution changed parameters
fluid(i).cnst.phi_p = 0.0300;
fluid(i).cnst.a_p = 8e-9;
fluid(i).cnst.mu_s = 0.41;
fluid(i).cnst.G_0 = 450;
fluid(i).cnst.sigma_y0 = 11;

fluid(i).par.W = 14.1442;
fluid(i).par.alfa = 0.4096;
fluid(i).par.b_0 = 0.0015;
fluid(i).par.d_f = 2.11;
fluid(i).par.porosity = 0.8192;
fluid(i).par.m_p = 468;

obj = fluid(i);
% Objective function
func = @(x) objectiveFunctionSilica(obj, x);

% Genetic algorithm
    options = optimoptions('ga','UseParallel',true, 'Display', 'iter');
    [x,fval,exitflag] = ga(func,numParameters,[],[],[],[],lb,ub,[],options);

% % Parallel tempering
%    x = ParTemp(func,parVec,1e-2,lb,'MIN',ub,'MAX');

% % Simulated annealing
% options = optimoptions('simulannealbnd', 'PlotFcns',...
%           {@saplotbestx,@saplotbestf,@saplotx,@saplotf},...
%           'OutputFcn', @myoutSA);
%    [x,fval,exitflag,output]  = simulannealbnd(func, lb+rand*(ub-lb), lb, ub, options);
    
% Assigning values to the parameters
    fluid(i).par.W = exp(x(1))-1;
    fluid(i).par.alfa = x(2);
    fluid(i).par.b_0 = exp(x(3));
    fluid(i).par.porosity = x(4);

    N = i;
    W = fluid(i).par.W;
    alfa =  fluid(i).par.alfa;
    b_0 =  fluid(i).par.b_0;
    porosity =  fluid(i).par.porosity;
    
    error = fval;
    
    if i ==1
    % Writing values to file
    fprintf("Writing the values to file\n");
    TT = table(N, W, alfa, b_0, porosity, error);

    writetable(TT,...
    "SAOptimizedSilica.csv",...
    'WriteVariableNames',1)
    else
       T2 = readtable('SAOptimizedSilica');
       T3 = table(N, W, alfa, b_0, porosity, error);
       TT = [T2;T3];
     writetable(TT,...
     "SAOptimizedSilica.csv",...
    'WriteVariableNames',1)
    end

end
