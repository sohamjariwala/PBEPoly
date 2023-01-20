%% Steady state
addpath("../")

par = [0.0053222798922441,0.958558875624029,3.90451317218712e-05,...
    2.36857200892795,1.00313330646949,302.212108530304,...
    1078.92571493154,10.599495769811,0.0459955288367757,...
    0.0224280851561805,0.00586678127507595]; %# Best fit

%% Initial PBEPoly class and load parameters 

obj = PBEPoly;

obj.par.W = par(1);
obj.par.alfa = par(2);
obj.par.b_0 = par(3);
obj.par.d_f = par(4);
obj.par.porosity = par(5);
obj.par.m_p = par(6);
obj.par.p = -2;
obj.par.kh = 0; 

obj.cnst.G_0 = par(7);
obj.cnst.sigma_y_0 = par(8);
obj.cnst.a_p = 20e-9;
obj.cnst.phi_p = par(9);
obj.cnst.mu_s = par(10);


HippetalSteadyShear

HippetalStepTransient

HippetalErrorCalculation
