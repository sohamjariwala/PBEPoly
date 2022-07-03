function fObj = objectiveFunctionSilica(obj, parVec)
% Objective function calculator for fitting the model to steady state and
% transient experimental data.

%% Steady state
    %% Loading the parameters 
    obj.par.W = exp(parVec(1))-1;
    obj.par.alfa = parVec(2);
    obj.par.b_0 = exp(parVec(3));
    obj.par.porosity = parVec(4);

    %% Loading the steady state experimental data set into variables
    silica_SS = readmatrix('silica.ExpData/silica_SS.txt');
    shear_rate_SS = silica_SS(:,1);
    shear_stress_SS = silica_SS(:,2);

    stress = zeros(size(shear_rate_SS));
    %% Solving the steady shear equations and collecting the output
    for i = length(shear_rate_SS):-1:1
        if i == length(shear_rate_SS)
            out = obj.steadyShear(shear_rate_SS(i));
            stress(i) = out.stress;
            logintMu(i,:) = out.logintMu;
        else
            out = obj.steadyShear(shear_rate_SS(i), out.logintMu);
            stress(i) = out.stress;
            logintMu(i,:) = out.logintMu;
        end
    end
    
    SS_error = norm((stress-shear_stress_SS)./(shear_stress_SS))...
        /length(shear_stress_SS);


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

% Clear temporary variables
clear opts

%% Solution
% Step down
    initial.EXITFLAG = 1;
    initial.logintMu = interp1(shear_rate_SS, logintMu, 5);
    initial.gamma_e = obj.gamma_lin;
    
    SD1 = stepShear(obj, iSD1, fSD1, time, initial);
    SD2 = stepShear(obj, iSD2, fSD2, time, initial);
    SD3 = stepShear(obj, iSD3, fSD3, time, initial);
    SD4 = stepShear(obj, iSD4, fSD4, time, initial);

% Step up
    initial.EXITFLAG = 1;
    initial.logintMu = interp1(shear_rate_SS, logintMu, 0.1);
    initial.gamma_e = obj.gamma_lin;
    
    SU1 = stepShear(obj, iSU1, fSU1, time, initial); 
    SU2 = stepShear(obj, iSU2, fSU2, time, initial);
    SU3 = stepShear(obj, iSU3, fSU3, time, initial);
    SU4 = stepShear(obj, iSU4, fSU4, time, initial);
    SU5 = stepShear(obj, iSU5, fSU5, time, initial);

%% Error
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

end