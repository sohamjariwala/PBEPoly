addpath('../');
%-------------------------------------------------------------------------------
obj = PBEPoly;

par = [19.89042756	0.622792833	0.002121292	2.281027346	0.857813642	582.6077229 ...
    35.48701115	0.490541689	0.35440899	0.424236925];

%-------------------------------------------------------------------------------
%% Loading the parameters
obj.par.W = par(1);
obj.par.alfa = par(2);
obj.par.b_0 = par(3);
obj.par.d_f = par(4);
obj.par.porosity = par(5);
obj.par.m_p = par(6);
obj.cnst.G_0 = par(7);
obj.cnst.sigma_y0 = par(8);
obj.cnst.mu_s = par(9);
obj.par.kh = par(10);
obj.par.p = 3;

%-------------------------------------------------------------------------------
%% Performing calculations 
loadExperimentalData;

WeietalSteadyShear; % Plot Fit

WeietalStepDownTransient;  % Plot Fits

WeietalShearReversal;

writeDataToOrigin

%% Error calculation
WeietalErrorCalculation;