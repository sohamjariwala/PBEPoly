function fObj = objectiveFunctionWei(obj, parVec)
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
    obj.par.kh = parVec(10);

    obj.par
    obj.cnst
    
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

%% Transient step shear plot
initial.EXITFLAG = 1;
initial.logintMu = interp1(SSEXP.shear_rate, logintMu, 0.1);
initial.stress = interp1(SSEXP.shear_rate, stress, 0.1);
initial.A = 1;

%% Set Step Shear parameters
i1 = 0.1; f1 = 1;
i2 = 0.1; f2 = 2.5;
i3 = 0.1; f3 = 5;

%% Solution
i1f1 = stepShear(obj, i1, f1, StepDownEXP.i0p1f1_t, initial);
i2f2 = stepShear(obj, i2, f2, StepDownEXP.i0p1f2p5_t, initial);
i3f3 = stepShear(obj, i3, f3, StepDownEXP.i0p1f5_t, initial);

transient_error = (norm((i1f1.stress-StepDownEXP.i0p1f1_stress)./mean(StepDownEXP.i0p1f1_stress)) + ...
                   norm((i2f2.stress-StepDownEXP.i0p1f2p5_stress)./mean(StepDownEXP.i0p1f2p5_stress)) + ...
                   norm((i3f3.stress-StepDownEXP.i0p1f5_stress)./mean(StepDownEXP.i0p1f5_stress)))...
                 ./length(StepDownEXP.i0p1f5_stress)/3;
%% Error
fObj = SS_error + transient_error;
end