%% Steady state
addpath("../")

par = [0.00532509034911643,0.968681143348009,4.5388112332297e-05, ...
    2.36906631541588,1.0083722449054,315.905368078551, ...
    1070.08332857968,10.4231736015752,0.045072528438335,0.0224376912336284];

%% Initial PBEPoly class and load parameters 

obj = PBEPoly;

obj.par.W = par(1);
obj.par.alfa = par(2);
obj.par.b_0 = par(3);
obj.par.d_f = par(4);
obj.par.porosity = par(5);
obj.par.m_p = par(6);
obj.par.p = 3;
obj.par.kh = 0; 

obj.cnst.G_0 = par(7);
obj.cnst.sigma_y_0 = par(8);
obj.cnst.a_p = 10e-9;
obj.cnst.phi_p = par(9);
obj.cnst.mu_s = par(10);

LoadExperimentalData

HippetalSteadyShear

% HippetalStepTransient

HippetalErrorCalculation
