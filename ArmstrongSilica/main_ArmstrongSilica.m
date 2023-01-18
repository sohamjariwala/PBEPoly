%% Loading the parameters and experimental data set into variables
addpath('../');
fluid = PBEPoly;

% par = [0.010508704	0.849999156	0.006340611	2.119061369	0.957866855	520.1646278 0];
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

silica_SS = readmatrix('silica.ExpData/silica_SS.txt');
silica_elastic_SS = readmatrix('silica.ExpData/silica_elastic_SS.txt');
shear_rate = silica_SS(:,1);
shear_stress_SS = silica_SS(:,2);
shear_rate_elastic = silica_elastic_SS(:,1);
elastic_comp_SS = silica_elastic_SS(:,2);

%% Transient experimental data
silica_stepDown = readmatrix('silica.ExpData/silica_Stepdown_i5.txt');
silica_stepDownTime = silica_stepDown(:,1);
silica_stepDowni5f2p5 = silica_stepDown(:,2);
silica_stepDowni5f1p0 = silica_stepDown(:,3);
silica_stepDowni5f0p5 = silica_stepDown(:,4);
silica_stepDowni5f0p1 = silica_stepDown(:,5); 
time = silica_stepDownTime;

% Set Step Shear parameters
tEnd = 100; %s
iSD1 = 5; fSD1 = 2.5;
iSD2 = 5; fSD2 = 1.0;
iSD3 = 5; fSD3 = 0.5;
iSD4 = 5; fSD4 = 0.1;
% time = [0 tEnd];

iSU1 = 0.1; fSU1 = 5.0;
iSU2 = 0.1; fSU2 = 2.5;
iSU3 = 0.1; fSU3 = 1.0;
iSU4 = 0.1; fSU4 = 0.5;
iSU5 = 0.1; fSU5 = 0.25;

opts = spreadsheetImportOptions("NumVariables", 6);

% Specify sheet and range
opts.Sheet = "Trans. Step Up Step Down Dat";
opts.DataRange = "F4:K53";

% Specify column names and types
opts.VariableNames = ["silica_stepUpTime", "silica_stepUpi0p1f5", "silica_stepUpi0p1f2p5", "silica_stepUpi0p1f1", "silica_stepUpi0p1f0p5", "silica_stepUpi0p1f0p25"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];

% Import the data
EXP = readtable("./silica.ExpData/2p9volx_fumed_silica_data_5apr15.xlsx", opts, "UseExcel", false);
silica_stepUpTime = EXP.silica_stepUpTime;
silica_stepUpi0p1f5 = EXP.silica_stepUpi0p1f5;
silica_stepUpi0p1f2p5 = EXP.silica_stepUpi0p1f2p5;
silica_stepUpi0p1f1 = EXP.silica_stepUpi0p1f1;
silica_stepUpi0p1f0p5 = EXP.silica_stepUpi0p1f0p5;
silica_stepUpi0p1f0p25 = EXP.silica_stepUpi0p1f0p25;

%% Clear temporary variables
clear opts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
%% Solving the steady shear equations and collecting the output variables
stress = zeros(size(shear_rate));
logintMu = (zeros(length(shear_rate), 5));

for i = length(shear_rate):-1:1
    if i == length(shear_rate)
        out = obj.steadyShear(shear_rate(i));
        out.A = 1;
    else
        out = obj.steadyShear(shear_rate(i), out);
    end
        EXITFLAG(i) = out.EXITFLAG;
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
        elastic_comp(i) = out.sigma_y;
        gamma_dot_p(i) = obj.gamma_dot_p(out.stress,1,out.logintMu,shear_rate(i));
        phi_max(i) = obj.phi_max(out.logintMu);
        rel_time(i) = obj.tau(out.logintMu);
end
toc;

%% Solving transient step shear equations
initial.EXITFLAG = 1;
initial.logintMu = interp1(shear_rate, logintMu, 5);
initial.stress = interp1(shear_rate, stress,5);
initial.A = 1;
% Step down
tic; SD1 = stepShear(obj, iSD1, fSD1, time, initial); toc;
tic; SD2 = stepShear(obj, iSD2, fSD2, time, initial); toc;
tic; SD3 = stepShear(obj, iSD3, fSD3, time, initial); toc;
tic; SD4 = stepShear(obj, iSD4, fSD4, time, initial); toc;

SD1_phi_a = zeros(size(time));
SD2_phi_a = SD1_phi_a;
SD3_phi_a = SD1_phi_a;
SD4_phi_a = SD1_phi_a;

for ii = 1:length(time)
    SD1_phi_a(ii) = obj.phi_a(SD1.logintMu(ii,:));
end
for ii = 1:length(time)
    SD2_phi_a(ii) = obj.phi_a(SD2.logintMu(ii,:));
end
for ii = 1:length(time)
   SD3_phi_a(ii) = obj.phi_a(SD3.logintMu(ii,:));
end
for ii = 1:length(time)
   SD4_phi_a(ii) = obj.phi_a(SD4.logintMu(ii,:));
end
% Step up
initial.EXITFLAG = 1;
initial.logintMu = interp1(shear_rate, logintMu, 0.1);
initial.stress = interp1(shear_rate, stress,0.1);

% Step up
tic; SU1 = stepShear(obj, iSU1, fSU1, time, initial); toc;
tic; SU2 = stepShear(obj, iSU2, fSU2, time, initial); toc;
tic; SU3 = stepShear(obj, iSU3, fSU3, time, initial); toc;
tic; SU4 = stepShear(obj, iSU4, fSU4, time, initial); toc;
tic; SU5 = stepShear(obj, iSU5, fSU5, time, initial); toc;

SU1_phi_a = zeros(size(time));
SU2_phi_a = SU1_phi_a;
SU3_phi_a = SU1_phi_a;
SU4_phi_a = SU1_phi_a;
SU5_phi_a = SU1_phi_a;

for ii = 1:length(time)
   SU1_phi_a(ii) = obj.phi_a(SU1.logintMu(ii,:));
end
for ii = 1:length(time)
   SU2_phi_a(ii) = obj.phi_a(SU2.logintMu(ii,:));
end
for ii = 1:length(time)
   SU3_phi_a(ii) = obj.phi_a(SU3.logintMu(ii,:));
end
for ii = 1:length(time)
   SU4_phi_a(ii) = obj.phi_a(SU4.logintMu(ii,:));
end
for ii = 1:length(time)
   SU5_phi_a(ii) = obj.phi_a(SU5.logintMu(ii,:));
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating radius of gyration
c = obj.MOMIC(logintMu);
for i=1:length(shear_rate)
    phi_a(i) = obj.phi_a(logintMu(i,:));
    c(:,i) = obj.MOMIC(logintMu(i,:));
    Rg(i) =obj.fra_moment(1+1/obj.par.d_f, c(:,i))*obj.cnst.a_p*10^9*2;
end
R_gOa_p = Rg/(obj.cnst.a_p*10^9);
% 
SS_error = norm((stress-shear_stress_SS)./(shear_stress_SS))...
        /length(shear_stress_SS);
    
transient_error_SD = ((norm((SD1.stress-silica_stepDowni5f2p5))./mean(silica_stepDowni5f2p5)) + ...
    (norm((SD2.stress-silica_stepDowni5f1p0))./mean(silica_stepDowni5f1p0)) + ...
    (norm((SD3.stress-silica_stepDowni5f0p5))./mean(silica_stepDowni5f0p5)) + ...
    (norm((SD4.stress-silica_stepDowni5f0p1))./mean(silica_stepDowni5f0p1)))  ...
    ./length(silica_stepDownTime)/4;

transient_error_SU = ((norm((SU1.stress-silica_stepUpi0p1f5))./mean(silica_stepUpi0p1f5)) + ...
    (norm((SU2.stress-silica_stepUpi0p1f2p5))./mean(silica_stepUpi0p1f2p5)) + ...
    (norm((SU3.stress-silica_stepUpi0p1f1))./mean(silica_stepUpi0p1f1)) + ...
    (norm((SU4.stress-silica_stepUpi0p1f0p5))./mean(silica_stepUpi0p1f0p5)) +...
    (norm((SU5.stress-silica_stepUpi0p1f0p25))./mean(silica_stepUpi0p1f0p25))) ...
    ./length(silica_stepUpTime)/5;


fObj = SS_error + transient_error_SD + transient_error_SU;
    
fprintf("Error in steady state fit = %f\n", SS_error);
fprintf("Transient step down error = %f\n", transient_error_SD);
fprintf("Transient step up error = %f\n", transient_error_SU);
fprintf("Total error = %1.15f\n", fObj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot for steady state
figure
loglog(shear_rate, stress,'LineWidth',2); hold on
loglog(shear_rate, shear_stress_SS, 's','LineWidth',2);
loglog(shear_rate, elastic_comp,'LineWidth',2);
loglog(shear_rate_elastic, elastic_comp_SS, '^', 'LineWidth',2);
loglog(shear_rate, stress - elastic_comp','LineWidth',2)
xlabel('Shear rate ($\mathrm{s}^{-1}$)','Interpreter','latex','FontSize',18);
ylabel('Stress (Pa)','Interpreter','latex','FontSize',18);
axis([-inf inf 0.1 inf])
yyaxis right
plot(shear_rate, phi_a, '--','LineWidth',2);
ylabel ( '$\phi_a$','Interpreter','latex','FontSize',18);

legend('Total Shear Stress (Model)', 'Total Shear Stress (Exp.)', ...
    'Elastic Component (Model)', 'Elastic Component (Exp.)', ...
    'Viscous Component (Model)', ...
    'Volume fraction (6 moments approximation)', ...
    'Interpreter','latex','FontSize',10,'location','best');

figure
semilogx(shear_rate,R_gOa_p, 'LineWidth',2);
xlabel('Shear rate ($\mathrm{s}^{-1}$)','Interpreter','latex','FontSize',18);
ylabel('$R_a/a_p$', 'Interpreter','latex','FontSize',18);
axis([-inf inf -inf inf])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots for transients
figure
box on;
semilogx(time, SD1.stress,'k-', ...
    time, SD2.stress,'r--',...
    time, SD3.stress,'m-.',...
    time, SD4.stress,'b:',...
    silica_stepDownTime, silica_stepDowni5f2p5, 'k^',...
    silica_stepDownTime, silica_stepDowni5f1p0, 'ro',...
    silica_stepDownTime, silica_stepDowni5f0p5, 'mv',...
    silica_stepDownTime, silica_stepDowni5f0p1, 'bs',...
    'MarkerSize',6,'LineWidth',2)

set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Time (s)');
ylabel(' Stress (Pa)');
legend('2.5 s^{-1}','1 s^{-1}','0.5 s^{-1}','0.1 s^{-1}');
axis([-inf inf 0 25])
grid on;

figure
box on;
semilogx(time, SD1_phi_a,'k-', ...
        time, SD2_phi_a,'r--',...
        time, SD3_phi_a,'m-.',...
        time, SD4_phi_a,'b:',...
        'MarkerSize',6,'LineWidth',2)

    set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Time (s)');
ylabel('\phi_a');
legend('2.5 s^{-1}','1 s^{-1}','0.5 s^{-1}', '0.1 s^{-1}');
axis([0.001 tEnd -inf inf])
grid on;

%-------------------------------------------------------------------------------
figure
box on;
semilogx(time, SU1.stress,'k-', ...
    time, SU2.stress,'r--',...
    time, SU3.stress,'m-.',...
    time, SU4.stress,'b:',...
    time, SU5.stress,'g.-',...
    silica_stepUpTime, silica_stepUpi0p1f5, 'k^',...
    silica_stepUpTime, silica_stepUpi0p1f2p5, 'ro',...
    silica_stepUpTime, silica_stepUpi0p1f1, 'mv',...
    silica_stepUpTime, silica_stepUpi0p1f0p5, 'bs',...
    silica_stepUpTime, silica_stepUpi0p1f0p25, 'gp',...
    'MarkerSize',6,'LineWidth',2)

set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Time (s)');
ylabel(' Stress (Pa)');
legend('5 s^{-1}','2.5 s^{-1}','0.1 s^{-1}','0.5 s^{-1}','2.5 s^{-1}');
axis([-inf inf 0 inf])
grid on;

figure
box on;
semilogx(time, SU1_phi_a,'k-', ...
        time, SU2_phi_a,'r--',...
        time, SU3_phi_a,'m-.',...
        time, SU4_phi_a,'b:',...
        time, SU5_phi_a,'g.-',...
        'MarkerSize',6,'LineWidth',2)

    set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Time (s)');
ylabel('\phi_a');
legend('5 s^{-1}','2.5 s^{-1}','1.0 s^{-1}','0.5 s^{-1}','0.25 s^{-1}');
axis([-inf inf -inf inf])
grid on;