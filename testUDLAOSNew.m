clear; 
testSilicaNew
close all;

%% Set UDLAOS parameters
omega = 1;
gamma_0 = 5;

%% Solution
nCycles = 5;
N = 1000;
initial.EXITFLAG = 1;
initial.logintMu = interp1(shear_rate,logintMu,gamma_0*omega);
% initial.logintMu = [-6.76767707468760,7.68560437728107,16.3664261808579,25.7075055546743,35.4709026817575];
initial.gamma_e = obj.gamma_lin;

time = linspace(0,2*pi*nCycles/omega,N)';
tic;out = UDLAOS(obj, gamma_0, omega, time, initial);toc

%% Post processing
strain = @(t) gamma_0*sin(omega*t);
strain_rate = @(t) gamma_0*omega + gamma_0*omega*cos(omega*t);

oscillatoryStrain = strain(linspace(0,2*nCycles*pi/omega,N));
oscillatoryStrainRate = strain_rate(linspace(0,2*nCycles*pi/omega,N));
np = N - floor(N/nCycles);
% np = 1;

%% Experimental data
silica_UDLAOS = readmatrix('silica.Expdata/silica_UDLAOS(1,5).txt');
silica_gamma_Exp = silica_UDLAOS(:,2)- gamma_0*omega*silica_UDLAOS(:,1);
silica_gamma_dot_Exp = silica_UDLAOS(:,3);
silica_sigma_Exp = silica_UDLAOS(:,4);

%% Plotting the results
figure
plot(oscillatoryStrain(np:end), out.stress(np:end)); hold on
plot(silica_gamma_Exp, silica_sigma_Exp, '.');
ylabel('Oscillatory shear stress (Pa)')
xlabel('Oscillatory strain (-)')

figure
plot(oscillatoryStrainRate(np:end), out.stress(np:end)); hold on
plot(silica_gamma_dot_Exp, silica_sigma_Exp, '.');
ylabel('Oscillatory shear stress (Pa)')
xlabel('Oscillatory strain rate (s^{-1})')


for i = 1: length(out.logintMu(:,1))
    phi_a(i) = obj.phi_a(out.logintMu(i,:));
end

figure
plot(oscillatoryStrainRate(np:end), phi_a(np:end));
xlabel('Oscillatory strain rate (s^{-1})')
ylabel('\phi_a (-)')

yyaxis right
plot(oscillatoryStrainRate(np:end), out.gamma_e(np:end));
