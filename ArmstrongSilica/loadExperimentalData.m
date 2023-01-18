%% Steady shear experimental data
silica_SS = readmatrix('silica.ExpData/silica_SS.txt');
silica_elastic_SS = readmatrix('silica.ExpData/silica_elastic_SS.txt');
shear_rate = silica_SS(:,1);
shear_stress_SS = silica_SS(:,2);
shear_rate_elastic = silica_elastic_SS(:,1);
elastic_comp_SS = silica_elastic_SS(:,2);

%% Step shear transient experimental data
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

clear opts

%% UDLAOS experimental data
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
