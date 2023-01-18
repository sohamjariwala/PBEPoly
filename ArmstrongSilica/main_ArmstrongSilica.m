%% Loading the parameters and experimental data set into variables
addpath('../');
fluid = PBEPoly;

par = [0.010988954	0.808530132	0.005761337	2.095445731	0.952727554	369.6420431	0];

% Changing constants and parameters to ones obtained from monodisperse
% solution changed parameters
obj = fluid;

% obj.par.W = exp(parVec(1))-1;
% obj.par.alfa = parVec(2);
% obj.par.b_0 = exp(parVec(3));
% obj.par.d_f = parVec(4);
% obj.par.porosity = parVec(5);
% obj.par.m_p = exp(parVec(6));

obj.par.W = par(1);
obj.par.alfa = par(2);
obj.par.b_0 = par(3);
obj.par.d_f = par(4);
obj.par.porosity = par(5);
obj.par.m_p = par(6);
obj.par.kh = par(7);
obj.par.p = 3;
obj.cnst.G_0 = 560;

loadExperimentalData

ArmstrongSilicaSteadyShear

ArmstrongSilicaStepDownTransient

ArmstrongSilicaStepUpTransient

ArmstrongSilicaUDLAOS

ArmstrongSilicaErrorCalculation