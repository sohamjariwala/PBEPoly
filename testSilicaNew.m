%% Loading the parameters and experimental data set into variables
fluid = PBEPoly;
% Changing constants and parameters to ones obtained from monodisperse 
% solution changed parameters
fluid.cnst.phi_p = 0.0300;
fluid.cnst.a_p = 8e-9;
fluid.cnst.mu_s = 0.41;
fluid.cnst.G_0 = 560;
fluid.cnst.sigma_y0 = 11;

fluid.par.W = 0.0380;
fluid.par.alfa = 1.2351;
fluid.par.b_0 = 0.004;

fluid.par.d_f = 2.11;
fluid.par.porosity = 8.973473867506506e-01;
fluid.par.m_p = 468;

obj = fluid;

% parVec = [3.92457341257613,0.609254640317527,-7.35735021647671,2.14020896364067,0.928691119044724,6.13420563599903];
% parVec = [3.88300750109340,0.604629636881591,-7.02580227865398,2.12242581912801,0.916637139213415,6.13621932587257];
% parVec = [4.01028119220008	0.781266793906812	-7.07782803745344	0.892993749305267];
% parVec = [3.30452    0.62327   -5.9527 0.9466];
% parVec = [4.017575 0.8641106 -5.908597 0.9411];
% parVec = [4.045 0.91 -5.70588 0.95];
% parVec = [4.0142 0.98 -5.8524 0.94];
% parVec = [4.04 0.98 -5.79844 0.9455];
% parVec = [4.235945210849526   1.207534999424150  -5.379382339397875   0.950512366546402];
% parVec = [1.1662 1.51458 -5.3424 0.89566];
% parVec = [3.64230556388067,0.620057132113812,-6.89776680411588,0.859826375826576];
% parVec = [0.072035895235900   1.179883317834147  -6.885865258517703   0.848891407123182];
% parVec = [0.068810601295548   0.924523717200031  -7.136220414585996  0.851765876896902];
% parVec = [4.1263 0.7213 -6.9142 0.8999];
% parVec = parVecBest;
% parVec = [0.016293920207754   1.649975101327514  -6.863277109520129   0.857693893803876];
% parVec = [0.0183, 2.04, -6.3073, 0.8471];
% 
%     obj.par.W = (exp(parVec(1))-1);
%     obj.par.alfa = parVec(2); % scaled by 0.5
%     obj.par.b_0 = exp(parVec(3)); % scaled by 0.5
%     obj.par.porosity = parVec(4);


X = [47.69896134 1 0.85*0.000981204 2.102684994 0.933224554 330.5396371 0.143827914];
X = [4.0869896134 0.61360438 0.0082 2.1102684994 0.923224554 468];
X = [40.71217784 0.978138515 0.002044432 0.934015942];

    obj.par.W = X(1);
    obj.par.alfa = X(2);
    obj.par.b_0 = X(3);
    obj.par.porosity = X(4);

    
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
        out = obj.steadyShearODE(shear_rate(i));
        EXITFLAG(i) = out.EXITFLAG;
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
        relaxation_time(i) = (obj.eta(obj.phi_a(out.logintMu),logintMu))...
                        ./obj.G(out.logintMu);
    else
        out = obj.steadyShearODE(shear_rate(i), out);
        EXITFLAG(i) = out.EXITFLAG;
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
        relaxation_time(i) = (obj.eta(obj.phi_a(out.logintMu),logintMu))...
                        ./obj.G(out.logintMu);
    end
end
toc;

% %% Solving transient step shear equations
% initial.EXITFLAG = 1;
% initial.logintMu = interp1(shear_rate, logintMu, 5);
% initial.gamma_e = obj.gamma_e_max(initial.logintMu);
% initial.stress = interp1(shear_rate, stress,5);
% % Step down
% tic; SD1 = stepShear(obj, iSD1, fSD1, time, initial); toc;
% tic; SD2 = stepShear(obj, iSD2, fSD2, time, initial); toc;
% tic; SD3 = stepShear(obj, iSD3, fSD3, time, initial); toc;
% tic; SD4 = stepShear(obj, iSD4, fSD4, time, initial); toc;
% 
% SD1_phi_a = zeros(size(time));
% SD2_phi_a = SD1_phi_a;
% SD3_phi_a = SD1_phi_a;
% SD4_phi_a = SD1_phi_a;
% 
% for ii = 1:length(time)
%     SD1_phi_a(ii) = obj.phi_a(SD1.logintMu(ii,:));
% end
% for ii = 1:length(time)
%     SD2_phi_a(ii) = obj.phi_a(SD2.logintMu(ii,:));
% end
% for ii = 1:length(time)
%    SD3_phi_a(ii) = obj.phi_a(SD3.logintMu(ii,:));
% end
% for ii = 1:length(time)
%    SD4_phi_a(ii) = obj.phi_a(SD4.logintMu(ii,:));
% end
% % Step up
% initial.EXITFLAG = 1;
% initial.logintMu = interp1(shear_rate, logintMu, 0.1);
% initial.gamma_e = obj.gamma_e_max(initial.logintMu);
% initial.stress = interp1(shear_rate, stress,0.1);
% 
% % Step down
% tic; SU1 = stepShear(obj, iSU1, fSU1, time, initial); toc;
% tic; SU2 = stepShear(obj, iSU2, fSU2, time, initial); toc;
% tic; SU3 = stepShear(obj, iSU3, fSU3, time, initial); toc;
% tic; SU4 = stepShear(obj, iSU4, fSU4, time, initial); toc;
% tic; SU5 = stepShear(obj, iSU5, fSU5, time, initial); toc;
% 
% SU1_phi_a = zeros(size(time));
% SU2_phi_a = SU1_phi_a;
% SU3_phi_a = SU1_phi_a;
% SU4_phi_a = SU1_phi_a;
% SU5_phi_a = SU1_phi_a;
% 
% for ii = 1:length(time)
%    SU1_phi_a(ii) = obj.phi_a(SU1.logintMu(ii,:));
% end
% for ii = 1:length(time)
%    SU2_phi_a(ii) = obj.phi_a(SU2.logintMu(ii,:));
% end
% for ii = 1:length(time)
%    SU3_phi_a(ii) = obj.phi_a(SU3.logintMu(ii,:));
% end
% for ii = 1:length(time)
%    SU4_phi_a(ii) = obj.phi_a(SU4.logintMu(ii,:));
% end
% for ii = 1:length(time)
%    SU5_phi_a(ii) = obj.phi_a(SU5.logintMu(ii,:));
% end
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating radius of gyration
c = obj.MOMIC(logintMu);
for i=1:length(shear_rate)
    phi_a(i) = obj.phi_a(logintMu(i,:));
    c(:,i) = obj.MOMIC(logintMu(i,:));
    Rg(i) =obj.fra_moment(1+1/obj.par.d_f, c(:,i))*obj.cnst.a_p*10^9*2;
    elastic_comp(i) = obj.elastic_stress(logintMu(i,:), obj.gamma_lin);
end
R_gOa_p = Rg/(obj.cnst.a_p*10^9);

SS_error = norm((stress-shear_stress_SS)./(shear_stress_SS))...
        /length(shear_stress_SS);
%     
% transient_error_SD = ((norm((SD1.stress-silica_stepDowni5f2p5))./mean(silica_stepDowni5f2p5)) + ...
%     (norm((SD2.stress-silica_stepDowni5f1p0))./mean(silica_stepDowni5f1p0)) + ...
%     (norm((SD3.stress-silica_stepDowni5f0p5))./mean(silica_stepDowni5f0p5)) + ...
%     (norm((SD4.stress-silica_stepDowni5f0p1))./mean(silica_stepDowni5f0p1)))  ...
%     ./length(silica_stepDownTime)/4;
% 
% transient_error_SU = ((norm((SU1.stress-silica_stepUpi0p1f5))./mean(silica_stepUpi0p1f5)) + ...
%     (norm((SU2.stress-silica_stepUpi0p1f2p5))./mean(silica_stepUpi0p1f2p5)) + ...
%     (norm((SU3.stress-silica_stepUpi0p1f1))./mean(silica_stepUpi0p1f1)) + ...
%     (norm((SU4.stress-silica_stepUpi0p1f0p5))./mean(silica_stepUpi0p1f0p5)) +...
%     (norm((SU5.stress-silica_stepUpi0p1f0p25))./mean(silica_stepUpi0p1f0p25))) ...
%     ./length(silica_stepUpTime)/5;
% 
% 
% fObj = SS_error + transient_error_SD + transient_error_SU;
%     
% fprintf("Error in steady state fit = %f\n", SS_error);
% fprintf("Transient step down error = %f\n", transient_error_SD);
% fprintf("Transient step up error = %f\n", transient_error_SU);
% fprintf("Total error = %f\n", fObj);
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
axis([0.01 inf -inf inf])
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
axis([0.01 inf -inf inf])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots for transients
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
semilogx(time, SD1.gamma_e,'k', ...
    time, SD2.gamma_e,'r',...
    time, SD3.gamma_e,'m',...
    time, SD4.gamma_e,'b',...
    'LineWidth',2);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Time (s)');
ylabel(' Elastic shear strain (-)');
legend('2.5 s^{-1}','1 s^{-1}','0.5 s^{-1}','0.1 s^{-1}');
axis([1e-2 1e3 0 inf])
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