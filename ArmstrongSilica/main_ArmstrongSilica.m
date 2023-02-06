%% Loading the parameters and experimental data set into variables
addpath('../');
fluid = PBEPoly;

par = [0.00910811	0.761403955	0.006881736	2.04518119	0.938844713	220.2195182	8.575409274	808.0949588];

% Changing constants and parameters to ones obtained from monodisperse
% solution changed parameters
obj = fluid;

obj.par.W = par(1);
obj.par.alfa = par(2);
obj.par.b_0 = par(3);
obj.par.d_f = par(4);
obj.par.porosity = par(5);
obj.par.m_p = par(6);
obj.par.kh = par(7);
obj.par.p = 3;
obj.cnst.G_0 = par(8);
obj.cnst.sigma_y0 = obj.cnst.sigma_y0 - obj.par.kh;

loadExperimentalData

ArmstrongSilicaSteadyShear

ArmstrongSilicaStepDownTransient

ArmstrongSilicaStepUpTransient

ArmstrongSilicaUDLAOS

ArmstrongSilicaLAOS

writeDataToOrigin

ArmstrongSilicaErrorCalculation