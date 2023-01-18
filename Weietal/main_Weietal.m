addpath('../');
%-------------------------------------------------------------------------------
obj = PBEPoly;

par = [20.99101415	0.606593006	0.00223603	2.22787788	0.858138915	336.6150101 ...
    36.73958518	0.795873736	0.338754722];

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
obj.par.kh = 0.1;
obj.par.p = 4;

%-------------------------------------------------------------------------------
%% Performing calculations 
loadExperimentalData;

WeietalSteadyShear; % Plot Fit

WeietalStepDownTransient;  % Plot Fits

% writeDataToOrigin

%% Error calculation
WeietalErrorCalculation;