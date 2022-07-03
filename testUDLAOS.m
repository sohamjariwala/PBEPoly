clear; close all;
%% Load parameters
load('silicaTest4.mat'); obj = silica;
obj = obj.inverseOfA;

% %     obj.par.W = 50; 
% %     obj.par.alfa = 0.61;
% %     obj.par.b_0 = 0.22e-2;
% 
%     obj.par.W = 0.008998930690468; 
%     obj.par.alfa = 0.665487794544652;
%     obj.par.b_0 = 0.003742494709682;

%% Set UDLAOS parameters
omega = 1;
gamma_0 = 5;

%% Experimental data
silica_UDLAOS = readmatrix('silica.Expdata/silica_UDLAOS(1,5).txt');
silica_gamma_Exp = silica_UDLAOS(:,2)- gamma_0*omega*silica_UDLAOS(:,1);
silica_gamma_dot_Exp = silica_UDLAOS(:,3);
silica_sigma_Exp = silica_UDLAOS(:,4);

%% Solution
nCycles = 5;
N = 1000;
tic; out = obj.UDLAOS(gamma_0, omega, linspace(0,2*nCycles*pi/omega,N)); toc;

%% Post processing
strain = @(t) gamma_0*sin(omega*t);
strain_rate = @(t) gamma_0*omega + gamma_0*omega*cos(omega*t);

oscillatoryStrain = strain(linspace(0,2*nCycles*pi/omega,N));
oscillatoryStrainRate = strain_rate(linspace(0,2*nCycles*pi/omega,N));
np = N - floor(N/nCycles);
% np = 1;
%% Plotting the results
figure(1)
plot(oscillatoryStrain(np:end), out.stress(np:end)); hold on
plot(silica_gamma_Exp, silica_sigma_Exp, '.');
ylabel('Oscillatory shear stress (Pa)')
xlabel('Oscillatory strain (-)')

figure(2)
plot(oscillatoryStrainRate(np:end), out.stress(np:end)); hold on
plot(silica_gamma_dot_Exp, silica_sigma_Exp, '.');
ylabel('Oscillatory shear stress (Pa)')
xlabel('Oscillatory strain rate (s^{-1})')


for i = 1: length(out.logintMu(:,1))
    phi_a(i) = obj.phi_a(out.logintMu(i,:));
end

figure(3)
plot(oscillatoryStrainRate(np:end), phi_a(np:end));
xlabel('Oscillatory strain rate (s^{-1})')
ylabel('\phi_a (-)')

yyaxis right
plot(oscillatoryStrainRate(np:end), out.gamma_e(np:end));
