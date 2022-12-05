parVec = [3.013704391154364   0.590198260682092  -6.354991898429272   2.268634081422561   0.883034591106096   5.953417691287351  37.277275082280944   0.782275208105718 0.333520407166337];
obj = PBEPoly;

% Changing constants and parameters to ones obtained from monodisperse 
% solution changed parameters
obj.cnst.phi_p = 0.0300;
obj.cnst.a_p = 8e-9;
obj.cnst.mu_s = 0.3592;
obj.cnst.G_0 = 2.3896;
obj.cnst.sigma_y0 = 0.7951;

obj.par.W = 14.1442;
obj.par.alfa = 0.4096;
obj.par.b_0 = 0.0015;
obj.par.d_f = 2.3325;
obj.par.porosity = 0.8192;
obj.par.m_p = 978;

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
    obj.cnst.mu_s = parVec(9);
    
%% Setup the Import Options and import the data for transients
opts = spreadsheetImportOptions("NumVariables", 6);

% Specify sheet and range
opts.Sheet = "Fig. 3b step rate tests";
opts.DataRange = "A4:F83";

% Specify column names and types
opts.VariableNames = ...
    ["i0p1f1_t", "i0p1f1_stress", ...
    "i0p1f2p5_t", "i0p1f2p5_stress",...
    "i0p1f5_t", "i0p1f5_stress"];
opts.VariableTypes = ...
    ["double", "double", "double", "double", "double", "double"];

% Import the data
StepDownEXP = readtable("./experimental_data.xlsx",opts, ...
                        "UseExcel", false);
clear opts

%% Setup the Import Options and import the data steady state
opts = spreadsheetImportOptions("NumVariables", 2);

% Specify sheet and range
opts.Sheet = "Fig. 3a steady state flow curve";
opts.DataRange = "A5:B32";

% Specify column names and types
opts.VariableNames = ["shear_rate", "stress"];
opts.VariableTypes = ["double", "double"];

% Import the data
SSEXP = readtable("./experimental_data.xlsx", opts,...
                             "UseExcel", false);

clear opts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shear_rate = SSEXP.shear_rate;
stress = zeros(size(SSEXP.shear_rate));
logintMu = (zeros(length(SSEXP.shear_rate), 5));

%% Solving the steady shear equations and collecting the output
for i = length(SSEXP.shear_rate):-1:1
    if i == length(SSEXP.shear_rate)
        out = obj.steadyShearODE(SSEXP.shear_rate(i));
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
        tau(i) = obj.tau(out.logintMu,SSEXP.shear_rate(i),obj.gamma_lin);
    else
        init = out;
        out = obj.steadyShearODE(SSEXP.shear_rate(i), init);
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
        tau(i) = obj.tau(out.logintMu,SSEXP.shear_rate(i),obj.gamma_lin);
    end
end

%% Relaxation time
figure
loglog(SSEXP.shear_rate, tau,LineWidth=2)
xlabel('Shear rate (s^{-1})')
ylabel('Relaxation time (s)')
grid on
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');

SS_error = norm((stress-SSEXP.stress)./mean(SSEXP.stress))...
    /length(SSEXP.stress);

fprintf("Steady state error = %f\n", SS_error);

%% Steady shear plot
figure
loglog(SSEXP.shear_rate, stress, SSEXP.shear_rate, SSEXP.stress,'o',...
    'MarkerSize',6,'LineWidth',2)
xlabel('Shear rate (s^{-1})','FontSize',18);
ylabel('Stress (Pa)','FontSize',18); 
grid on;       
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
axis([-inf inf -inf 100]);

distribution