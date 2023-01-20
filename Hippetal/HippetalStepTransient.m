%% Transients
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

%% Solution
i1f1 = stressResponse(obj, i1, vulcan_stepDown1.shear_rate, vulcan_stepDown1.time, initial);
i2f2 = stressResponse(obj, i2, vulcan_stepDown2.shear_rate, vulcan_stepDown2.time, initial);
i3f3 = stressResponse(obj, i3, vulcan_stepDown3.shear_rate, vulcan_stepDown3.time, initial);
i4f4 = stressResponse(obj, i4, vulcan_stepDown4.shear_rate, vulcan_stepDown4.time, initial);

%% Plots
figure
semilogx(vulcan_stepDown1.time, vulcan_stepDown1.stress,'s',vulcan_stepDown1.time, i1f1.stress); hold on;
semilogx(vulcan_stepDown2.time, vulcan_stepDown2.stress,'s',vulcan_stepDown2.time, i2f2.stress);
semilogx(vulcan_stepDown3.time, vulcan_stepDown3.stress,'s',vulcan_stepDown3.time, i3f3.stress);
semilogx(vulcan_stepDown4.time, vulcan_stepDown4.stress,'s',vulcan_stepDown4.time, i4f4.stress);