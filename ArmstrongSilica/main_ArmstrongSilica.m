%% Loading the parameters and experimental data set into variables
addpath('../');
fluid = PBEPoly;

par = [0.010790136	0.815485863	0.005700755	2.087541192	0.959594226	265.7177498	10.7115235];
par = [0.009554792	0.828975549	0.004909345	2.096963699	0.955374652	299.7560198	7.7];
par = [0.009214781	0.830996836	0.006742349	2.074355917	0.94418072	334.6654253	8.563311753 560];
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