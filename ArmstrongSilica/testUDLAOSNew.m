main_ArmstrongSilica
    
%% Transient experimental data
UDLAOS1_Exp = readmatrix('silica.Expdata/silica_UDLAOS(1,1).txt');
UDLAOS2_Exp = readmatrix('silica.Expdata/silica_UDLAOS(1,5).txt');
UDLAOS3_Exp = readmatrix('silica.Expdata/silica_UDLAOS(1,10).txt');

% Extract time
Exp_time1 = UDLAOS1_Exp(:,1);
Exp_time2 = UDLAOS2_Exp(:,1);
Exp_time3 = UDLAOS3_Exp(:,1);

% Extract stress
Exp_stress1 = UDLAOS1_Exp(:,end);
Exp_stress2 = UDLAOS2_Exp(:,end);
Exp_stress3 = UDLAOS3_Exp(:,end);

omega1 = 1; gamma_01 = 1;
omega2 = 1; gamma_02 = 5;
omega3 = 1; gamma_03 = 10;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SS_error = norm((stress-shear_stress_SS)./(shear_stress_SS))...
        /length(shear_stress_SS);
    
transient_error_UDLAOS = ...
    ((norm((UDLAOS1.stress-Exp_stress1)./mean(Exp_stress1)))./length(Exp_stress1) + ...
    (norm((UDLAOS2.stress-Exp_stress2)./mean(Exp_stress2)))./length(Exp_stress2) + ...
    (norm((UDLAOS3.stress-Exp_stress3)./mean(Exp_stress3)))./length(Exp_stress3)) ...
    /3;

fObj = SS_error + transient_error_UDLAOS

%% Plot
% reduced strain

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
