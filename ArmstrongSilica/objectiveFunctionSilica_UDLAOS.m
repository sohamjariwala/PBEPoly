function fObj = objectiveFunctionSilica_UDLAOS(obj, parVec)
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

%% Solving transient UDLAOS equations
initial1.EXITFLAG = 1;
initial1.logintMu = interp1(shear_rate, logintMu, gamma_01*omega1);
initial1.stress = interp1(shear_rate, stress,gamma_01*omega1);
initial1.A = 1;

initial2.EXITFLAG = 1;
initial2.logintMu = interp1(shear_rate, logintMu, gamma_02*omega2);
initial2.stress = interp1(shear_rate, stress,gamma_02*omega2);
initial2.A = 1;

initial3.EXITFLAG = 1;
initial3.logintMu = interp1(shear_rate, logintMu, gamma_03*omega3);
initial3.stress = interp1(shear_rate, stress,gamma_03*omega3);
initial3.A = 1;


% UDLAOS
UDLAOS1 = UDLAOS(obj, gamma_01, omega1, Exp_time1, initial1);
UDLAOS2 = UDLAOS(obj, gamma_02, omega2, Exp_time2, initial2);
UDLAOS3 = UDLAOS(obj, gamma_03, omega3, Exp_time3, initial3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SS_error = norm((stress-shear_stress_SS)./(shear_stress_SS))...
        /length(shear_stress_SS);
    
transient_error_UDLAOS = ...
    ((norm((UDLAOS1.stress-Exp_stress1)./mean(Exp_stress1)))./length(Exp_stress1) + ...
    (norm((UDLAOS2.stress-Exp_stress2)./mean(Exp_stress2)))./length(Exp_stress2) + ...
    (norm((UDLAOS3.stress-Exp_stress3)./mean(Exp_stress3)))./length(Exp_stress3)) ...
    /3;

fObj = SS_error + transient_error_UDLAOS;

end