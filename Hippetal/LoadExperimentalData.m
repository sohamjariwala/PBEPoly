%% Loading the steady state experimental data set into variables
SSEXP.shear_rate = readmatrix('vulcan.ExpData/shear_rate_SS.txt');
SSEXP.stress = readmatrix('vulcan.ExpData/shear_stress_SS.txt');

RGEXP.shear_rate = readmatrix('vulcan.ExpData/shear_rate_RG.txt');
RGEXP.RG = readmatrix('vulcan.ExpData/RG.txt');

%% Experimental data
vulcan_stepDown1.stress = readmatrix('vulcan.ExpData/stress_StepDown_i2500f50.txt');
vulcan_stepDown2.stress = readmatrix('vulcan.ExpData/stress_StepDown_i2500f100.txt');
vulcan_stepDown3.stress = readmatrix('vulcan.ExpData/stress_StepDown_i2500f500.txt');
vulcan_stepDown4.stress = readmatrix('vulcan.ExpData/stress_StepDown_i2500f1000.txt');

vulcan_stepDown1.time = readmatrix('vulcan.ExpData/t_StepDown_i2500f50.txt');
vulcan_stepDown2.time = readmatrix('vulcan.ExpData/t_StepDown_i2500f100.txt');
vulcan_stepDown3.time = readmatrix('vulcan.ExpData/t_StepDown_i2500f500.txt');
vulcan_stepDown4.time = readmatrix('vulcan.ExpData/t_StepDown_i2500f1000.txt');

vulcan_stepDown1.shear_rate = readmatrix('vulcan.ExpData/shearRate_StepDown_i2500f50.txt');
vulcan_stepDown2.shear_rate = readmatrix('vulcan.ExpData/shearRate_StepDown_i2500f100.txt');
vulcan_stepDown3.shear_rate = readmatrix('vulcan.ExpData/shearRate_StepDown_i2500f500.txt');
vulcan_stepDown4.shear_rate = readmatrix('vulcan.ExpData/shearRate_StepDown_i2500f1000.txt');
