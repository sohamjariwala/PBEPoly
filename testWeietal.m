% parVec = [2.9378 0.570755 -6.0671 2.1708 0.90804 6.065031 39.067041 0.7951 0.1963722];
% parVec = [2.9008 0.5717 -5.9283 2.1635 0.9076 6.0341 38.2842 0.7951 0.2010];
% parVec = [2.864550050745758   0.572115011942363  -5.941326101945536   2.153231558815501   0.902963313001152   5.966540711533419  39.463950524306519   0.776901649455502   0.200855597641966];
% parVec = [2.803492926108994   0.572657166195471  -6.016251369726170   2.274097339535812   0.925136106766159   6.099712107958716  39.305903964253972   0.88217387338744   0.201662260698935];
% parVec = [3.0091 0.5748 -6.1299 2.28465 0.8962 5.9672 38.3639 0.7715 0.3481];
% parVec = [3.009184927963146   0.574820100362079  -6.129916683018121   2.284651220139327   0.896239432060398   5.967239794876706  38.363935423635844   0.771529572427426   0.348133957695984];
% parVec = [3.024952387572246   0.589041528128687  -6.249040808649854   2.270036781987128   0.887832090476037   5.840684821201803  38.501763148417659   0.790236543136271   0.339578696109709];
% parVec = [3.058294210350837   0.590840162774347  -6.279278825175675   2.265820000365196   0.883692869941822   5.888818412949577  38.135980833025492   0.790651168745462   0.336344273576225];
% parVec = [3.052969e+00, 5.902315e-01, -6.315614e+00, 2.265117e+00, 8.821384e-01, 5.878319e+00, 3.806119e+01, 7.888523e-01, 3.362045e-01];
% parVec = [3.018240e+00, 5.848876e-01, -6.412818e+00, 2.272982e+00, 8.816055e-01, 5.944847e+00, 3.712128e+01, 7.804288e-01, 3.360543e-01];
% parVec = [3.032708665273683   0.586249427294304  -6.387371315598383   2.273004409469752   0.882737572027301   5.969332515797490 37.111622267827748   0.782155476209003   0.334940333820390];
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
    
% Objective function calculator for fitting the model to steady state and
% transient experimental data.
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
        out = obj.steadyShear(SSEXP.shear_rate(i));
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
    else
        out = obj.steadyShear(SSEXP.shear_rate(i), out.logintMu);
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
    end
end

SS_error = norm((stress-SSEXP.stress)./mean(SSEXP.stress))...
    /length(SSEXP.stress);

fprintf("Steady state error = %f\n", SS_error);

%% Steady shear plot 
loglog(SSEXP.shear_rate, stress, SSEXP.shear_rate, SSEXP.stress,'o',...
    'MarkerSize',6,'LineWidth',2)
xlabel('Shear rate (s^{-1})','FontSize',18);
ylabel('Stress (Pa)','FontSize',18); 
grid on;       
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
axis([-inf inf 5e-1 100]);

distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transient step shear plot
initial.EXITFLAG = 1;
initial.logintMu = interp1(SSEXP.shear_rate, logintMu, 0.1);
initial.gamma_e = obj.gamma_lin;

    %% Set Step Shear parameters
    i1 = 0.1; f1 = 1;
    i2 = 0.1; f2 = 2.5;
    i3 = 0.1; f3 = 5;
    
    %% Solution
    i1f1 = stepShear(obj, i1, f1, StepDownEXP.i0p1f1_t, initial);
    i2f2 = stepShear(obj, i2, f2, StepDownEXP.i0p1f2p5_t, initial);
    i3f3 = stepShear(obj, i3, f3, StepDownEXP.i0p1f5_t, initial);

    transient_error = (norm((i1f1.stress-StepDownEXP.i0p1f1_stress)./(StepDownEXP.i0p1f1_stress)) + ...
                       norm((i2f2.stress-StepDownEXP.i0p1f2p5_stress)./(StepDownEXP.i0p1f2p5_stress)) + ...
                       norm((i3f3.stress-StepDownEXP.i0p1f5_stress)./(StepDownEXP.i0p1f5_stress)))...
                     ./length(StepDownEXP.i0p1f5_stress);

    fprintf("Transient step shear error = %f\n", transient_error);

 %% Transient plots
 figure
    semilogx(StepDownEXP.i0p1f1_t, StepDownEXP.i0p1f1_stress,'ko',...     
             StepDownEXP.i0p1f2p5_t, StepDownEXP.i0p1f2p5_stress, 'bv', ...
             StepDownEXP.i0p1f5_t, StepDownEXP.i0p1f5_stress, 'g^',...
             StepDownEXP.i0p1f5_t, i3f3.stress, 'g',...
             StepDownEXP.i0p1f1_t, i1f1.stress,'k',...
             StepDownEXP.i0p1f2p5_t, i2f2.stress, 'b', ...
            'MarkerSize',6,'LineWidth',2)
         
             legend('1 (1/s)','2.5 (1/s)','5 (1/s)');
    set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
    xlabel('Time (s)');
    ylabel(' Stress (Pa)');

    axis([-inf inf 0 12])
    grid on;       

%% Error
fObj = SS_error + transient_error;
fprintf("Total error = %f\n", fObj);