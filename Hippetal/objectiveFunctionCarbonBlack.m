function fObj = objectiveFunctionCarbonBlack(obj, parVec)
% Objective function calculator for fitting the model to steady state and
% transient experimental data.

%% Steady state
%% Loading the parameters 
    obj.par.W = exp(parVec(1))-1;
    obj.par.alfa = parVec(2);
    obj.par.b_0 = exp(parVec(3));
    obj.par.d_f = parVec(4);
    obj.par.porosity = parVec(5);
    obj.par.m_p = exp(parVec(6));
    obj.cnst.G_0 = parVec(7);
    obj.cnst.sigma_y0 = parVec(8);
    obj.cnst.phi_p = parVec(9);
    obj.cnst.mu_s = parVec(10);

%% Loading the steady state experimental data set into variables
SSEXP.shear_rate = readmatrix('vulcan.ExpData/shear_rate_SS.txt');
SSEXP.stress = readmatrix('vulcan.ExpData/shear_stress_SS.txt');

stress = zeros(size(SSEXP.shear_rate));
logintMu = (zeros(length(SSEXP.shear_rate), 5));

%% Solving the steady shear equations and collecting the output
for i = length(SSEXP.shear_rate):-1:1
    if i == length(SSEXP.shear_rate)
        out = obj.steadyShear(SSEXP.shear_rate(i));
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
    else
        out = obj.steadyShear(SSEXP.shear_rate(i), out);
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
    end
end

SS_error = norm((stress-SSEXP.stress)./(SSEXP.stress))...
    /length(SSEXP.stress);

%% Transient state
initial.EXITFLAG = 1;
initial.logintMu = interp1(SSEXP.shear_rate, logintMu, 2500);
initial.stress = interp1(SSEXP.shear_rate, stress, 2500);
initial.A = 1;

%% Loading the transient experimental data into variables
%% Set Step Shear parameters
i1 = 2500; f1 = 50;
i2 = 2500; f2 = 100;
i3 = 2500; f3 = 500;
i4 = 2500; f4 = 1000;

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

%% Solution
i1f1 = stepShear(obj, i1, f1, vulcan_stepDown1.time, initial);
i2f2 = stepShear(obj, i2, f2, vulcan_stepDown2.time, initial);
i3f3 = stepShear(obj, i3, f3, vulcan_stepDown3.time, initial);
i4f4 = stepShear(obj, i4, f4, vulcan_stepDown4.time, initial);

transient_error = ((norm((i1f1.stress-vulcan_stepDown1.stress)./mean(vulcan_stepDown1.stress))./length(vulcan_stepDown1.stress)) + ...
    (norm((i2f2.stress-vulcan_stepDown2.stress)./mean(vulcan_stepDown2.stress))./length(vulcan_stepDown2.stress)) + ...                   sum((norm(i3f3.stress-vulcan_stepDown3.stress))./mean(vulcan_stepDown3.stress)) + ...
    (norm((i3f3.stress-vulcan_stepDown3.stress)./mean(vulcan_stepDown3.stress))./length(vulcan_stepDown3.stress)) + ...
    (norm((i4f4.stress-vulcan_stepDown4.stress)./mean(vulcan_stepDown4.stress))./length(vulcan_stepDown4.stress)))/4;

%% Error
fObj = SS_error;
fObj = SS_error + transient_error;
end