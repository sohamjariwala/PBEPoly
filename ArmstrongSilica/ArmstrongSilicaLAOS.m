%% LAOS
nCycles = 12;

% omega = 0.01
initial1.EXITFLAG = 1;
initial1.logintMu = interp1(shear_rate, logintMu, 0.01);
initial1.stress = interp1(shear_rate, stress, 0.01);
initial1.A = 1;

time = linspace(0, nCycles*pi, 1000)/0.01;

LAOS0p01_1 = LAOS(obj, 1, 0.01, time, initial1);
LAOS0p01_10 = LAOS(obj, 10, 0.01, time, initial1);
LAOS0p01_100 = LAOS(obj, 100, 0.01, time, initial1);
LAOS0p01_1000 = LAOS(obj, 1000, 0.01, time, initial1);

% omega = 0.1
initial1.EXITFLAG = 1;
initial1.logintMu = interp1(shear_rate, logintMu, 0.1);
initial1.stress = interp1(shear_rate, stress,0.1);
initial1.A = 1;

time = linspace(0, nCycles*pi, 1000)/0.1;

LAOS0p1_0p1 = LAOS(obj, 0.1, 0.1, time, initial1);
LAOS0p1_1 = LAOS(obj, 1, 0.1, time, initial1);
LAOS0p1_10 = LAOS(obj, 10, 0.1, time, initial1);
LAOS0p1_100 = LAOS(obj, 100, 0.1, time, initial1);
LAOS0p1_1000 = LAOS(obj, 1000, 0.1, time, initial1);

% omega = 1
initial1.EXITFLAG = 1;
initial1.logintMu = interp1(shear_rate, logintMu, 1);
initial1.stress = interp1(shear_rate, stress, 1);
initial1.A = 1;

time = linspace(0, nCycles*pi, 1000)/1;

LAOS1_0p1 = LAOS(obj, 0.1, 1, time, initial1);
LAOS1_1 = LAOS(obj, 1, 1, time, initial1);
LAOS1_10 = LAOS(obj, 10, 1, time, initial1);
LAOS1_100 = LAOS(obj, 100, 1, time, initial1);

% omega = 10
initial1.EXITFLAG = 1;
initial1.logintMu = interp1(shear_rate, logintMu, 10);
initial1.stress = interp1(shear_rate, stress,10);
initial1.A = 1;

time = linspace(0, nCycles*pi, 1000)/10;

LAOS10_0p1 = LAOS(obj, 0.1, 10, time, initial1);
LAOS10_1 = LAOS(obj, 1, 10, time, initial1);
LAOS10_10 = LAOS(obj, 10, 10, time, initial1);