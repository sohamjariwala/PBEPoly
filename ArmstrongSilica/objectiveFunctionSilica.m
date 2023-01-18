function fObj = objectiveFunctionSilica(obj, parVec)
% Objective function calculator for fitting the model to steady state and
% transient experimental data.

%% Steady state
%% Loading the parameters 
obj.par.W = exp(parVec(1))-1;
obj.par.alfa = parVec(2);
obj.par.b_0 = exp(parVec(3));
obj.par.d_f = parVec(4);
obj.par.porosity = parVec(5);
obj.par.m_p = exp(parVec(6));
obj.par.kh = parVec(7);

%% Steady state experimental data    
silica_SS = readmatrix('silica.ExpData/silica_SS.txt');
shear_rate = silica_SS(:,1);
shear_stress_SS = silica_SS(:,2);
    
%% Transient experimental data
silica_stepDown = readmatrix('silica.ExpData/silica_Stepdown_i5.txt');
silica_stepDownTime = silica_stepDown(:,1);
silica_stepDowni5f2p5 = silica_stepDown(:,2);
silica_stepDowni5f1p0 = silica_stepDown(:,3);
silica_stepDowni5f0p5 = silica_stepDown(:,4);
silica_stepDowni5f0p1 = silica_stepDown(:,5); 
time = silica_stepDownTime;

% Set Step Shear parameters
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

%% UDLAOS
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
end


%% Solving transient step shear equations
initial.EXITFLAG = 1;
initial.logintMu = interp1(shear_rate, logintMu, 5);
initial.stress = interp1(shear_rate, stress,5);
initial.A = 1;

% Step down
SD1 = stepShear(obj, iSD1, fSD1, time, initial);
SD2 = stepShear(obj, iSD2, fSD2, time, initial);
SD3 = stepShear(obj, iSD3, fSD3, time, initial);
SD4 = stepShear(obj, iSD4, fSD4, time, initial);

% Step up
initial.EXITFLAG = 1;
initial.logintMu = interp1(shear_rate, logintMu, 0.1);
initial.stress = interp1(shear_rate, stress,0.1);

% Step up
SU1 = stepShear(obj, iSU1, fSU1, time, initial);
SU2 = stepShear(obj, iSU2, fSU2, time, initial);
SU3 = stepShear(obj, iSU3, fSU3, time, initial);
SU4 = stepShear(obj, iSU4, fSU4, time, initial);
SU5 = stepShear(obj, iSU5, fSU5, time, initial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Solving transient UDLAOS equations
% initial1.EXITFLAG = 1;
% initial1.logintMu = interp1(shear_rate, logintMu, gamma_01*omega1);
% initial1.stress = interp1(shear_rate, stress,gamma_01*omega1);
% initial1.A = 1;
% 
% initial2.EXITFLAG = 1;
% initial2.logintMu = interp1(shear_rate, logintMu, gamma_02*omega2);
% initial2.stress = interp1(shear_rate, stress,gamma_02*omega2);
% initial2.A = 1;
% 
% initial3.EXITFLAG = 1;
% initial3.logintMu = interp1(shear_rate, logintMu, gamma_03*omega3);
% initial3.stress = interp1(shear_rate, stress,gamma_03*omega3);
% initial3.A = 1;
% 
% 
% % UDLAOS
% UDLAOS1 = UDLAOS(obj, gamma_01, omega1, Exp_time1, initial1);
% UDLAOS2 = UDLAOS(obj, gamma_02, omega2, Exp_time2, initial2);
% UDLAOS3 = UDLAOS(obj, gamma_03, omega3, Exp_time3, initial3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% transient_error_UDLAOS = ...
%     ((norm((UDLAOS1.stress-Exp_stress1)./mean(Exp_stress1)))./length(Exp_stress1) + ...
%     (norm((UDLAOS2.stress-Exp_stress2)./mean(Exp_stress2)))./length(Exp_stress2) + ...
%     (norm((UDLAOS3.stress-Exp_stress3)./mean(Exp_stress3)))./length(Exp_stress3)) ...
%     /3;

fObj = SS_error + transient_error_SD + transient_error_SU;

end