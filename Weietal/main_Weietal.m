addpath('../');
%-------------------------------------------------------------------------------
parVec = [3.091690100862115   0.604932033000000  -6.177390351405440 ...
    2.300037015000000   0.879680247000000   5.732836850990885 ...
    36.539401529999999   0.795985449000000   0.323856153000000];
obj = PBEPoly;

par = [21.68057883	0.584237984	0.003467759	2.24475337	0.886963815	298.725812 ...
    35.89360516	0.859146331	0.304398054];


%% Steady state
% Loading the parameters 
obj.par.W = par(1);
obj.par.alfa = par(2);
obj.par.b_0 = par(3);
obj.par.d_f = par(4);
obj.par.porosity = par(5);
obj.par.m_p = par(6);
obj.cnst.G_0 = par(7);
obj.cnst.sigma_y0 = par(8);
obj.cnst.mu_s = par(9);
obj.par.p = 4;

%-------------------------------------------------------------------------------
%% Performing calculations 
loadExperimentalData;

WeietalSteadyShear; % Fit

WeietalStepDownTransient;  % Fit

% writeDataToOrigin

%% Error
fObj = SS_error + transient_error;
fprintf("Total error = %f\n", fObj);