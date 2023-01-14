% Objective function calculator for fitting the model to steady state and
% transient experimental data.

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

%% Setup the Import Options and import the data for Step down transients
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

%% Step the Import Options and import the data for Step Down transients
opts = spreadsheetImportOptions("NumVariables", 2);

opts.Sheet = "Fig. 4";
opts.DataRange = "A5:B59";

opts.VariableNames = ...
    ["i3f9_t","i3f9_stress"];
opts.VariableTypes = ["double","double"];

StepDownEXP_Test = readtable("./experimental_data.xlsx",opts,"UseExcel",false);




