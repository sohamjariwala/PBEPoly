addpath('../');
%-------------------------------------------------------------------------------
parVec = [3.091690100862115   0.604932033000000  -6.177390351405440 ...
    2.300037015000000   0.879680247000000   5.732836850990885 ...
    36.539401529999999   0.795985449000000   0.323856153000000];
obj = PBEPoly;

par = [21.68057883	0.584237984	0.003467759	2.24475337	0.886963815	298.725812 ...
    35.89360516	0.859146331	0.304398054];

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

writeDataToOrigin

%% Error calculation
% WeietalErrorCalculation;