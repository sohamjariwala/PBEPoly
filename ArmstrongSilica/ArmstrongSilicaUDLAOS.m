%% Solving transient UDLAOS equations
initial1.EXITFLAG = 1;
initial1.logintMu = interp1(shear_rate, logintMu, gamma_01*omega1);
initial1.stress = interp1(shear_rate, stress,gamma_01*omega1);
initial1.A = 1;

initial2.EXITFLAG = 1;
initial2.logintMu = interp1(shear_rate, logintMu, gamma_02*omega2);
initial2.stress = interp1(shear_rate, stress,gamma_02*omega2);
initial2.A = 1;

initial3.EXITFLAG = 1;
initial3.logintMu = interp1(shear_rate, logintMu, gamma_03*omega3);
initial3.stress = interp1(shear_rate, stress,gamma_03*omega3);
initial3.A = 1;

% UDLAOS
UDLAOS1 = UDLAOS(obj, gamma_01, omega1, Exp_time1, initial1);
UDLAOS2 = UDLAOS(obj, gamma_02, omega2, Exp_time2, initial2);
UDLAOS3 = UDLAOS(obj, gamma_03, omega3, Exp_time3, initial3);
    
%% Plot

figure('Name', 'UDLAOS | Elastic projection')
reducedStrain = @(t, gamma_0, omega) gamma_0*sin(omega*t);

plot(reducedStrain(Exp_time1, gamma_01, omega1), Exp_stress1, 'k', ...
    reducedStrain(Exp_time1, gamma_01, omega1), UDLAOS1.stress,'r', ...
    'LineWidth',2); hold on;
plot(reducedStrain(Exp_time2, gamma_02, omega2), Exp_stress2, 'k', ...
    reducedStrain(Exp_time2, gamma_02, omega2), UDLAOS2.stress,'r', ...
    'LineWidth',2);
plot(reducedStrain(Exp_time3, gamma_03, omega3), Exp_stress3, 'k', ...
    reducedStrain(Exp_time3, gamma_03, omega3), UDLAOS3.stress,'r', ...
    'LineWidth',2);

xlabel('Oscillatory strain (-)');
ylabel('Stress (Pa)');
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');

%-------------------------------------------------------------------------------

figure('Name', 'UDLAOS | Viscous projection')
reducedStrainRate = @(t, gamma_0, omega) gamma_0*omega + gamma_0*omega*cos(omega*t);
plot(reducedStrainRate(Exp_time1, gamma_01, omega1), Exp_stress1, 'k', ...
    reducedStrainRate(Exp_time1, gamma_01, omega1), UDLAOS1.stress,'r', ...
    'LineWidth',2); hold on;
plot(reducedStrainRate(Exp_time2, gamma_02, omega2), Exp_stress2, 'k', ...
    reducedStrainRate(Exp_time2, gamma_02, omega2), UDLAOS2.stress,'r', ...
    'LineWidth',2);
plot(reducedStrainRate(Exp_time3, gamma_03, omega3), Exp_stress3, 'k', ...
    reducedStrainRate(Exp_time3, gamma_03, omega3), UDLAOS3.stress,'r', ...
    'LineWidth',2);
xlabel('Osciallatory strain rate (s^{-1})');
ylabel('Stress (Pa)');
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
