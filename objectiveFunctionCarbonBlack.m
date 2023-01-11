function fObj = objectiveFunctionCarbonBlack(obj, parVec)
% Objective function calculator for fitting the model to steady state and
% transient experimental data.

%% Steady state
    %% Loading the parameters 
    obj.par.W = exp(parVec(1))-1;
    obj.par.alfa = parVec(2);
    obj.par.b_0 = exp(parVec(3));

    %% Loading the steady state experimental data set into variables
    silica_SS = readmatrix('silica.ExpData/silica_SS.txt');
    shear_rate_SS = silica_SS(:,1);
    shear_stress_exp = silica_SS(:,2);

    stress = zeros(size(shear_rate_SS));
    %% Solving the steady shear equations and collecting the output
    for i = length(shear_rate_SS):-1:1
        if i == length(shear_rate_SS)
            out = obj.steadyShear(shear_rate_SS(i));
            SS_stress(i) = out.stress;
        else
            out = obj.steadyShear(shear_rate_SS(i), out.logintMu);
            SS_stress(i) = out.stress;
        end
    end

    SS_error = sum(norm(stress-shear_stress_exp)./mean(shear_stress_exp))...
        /length(shear_stress_exp);


%% Transient state
    %% Loading the transient experimental data into variables
    %% Set Step Shear parameters
    i1 = 5; f1 = 2.5;
    i2 = 5; f2 = 1.0;
    i3 = 5; f3 = 0.5;
    i4 = 5; f4 = 0.1;

    %% Experimental data
    silica_stepDown = readmatrix('silica.ExpData/silica_Stepdown_i5.txt');
    silica_stepDownTime = silica_stepDown(:,1);
    silica_stepDowni5f2p5 = silica_stepDown(:,2);
    silica_stepDowni5f1p0 = silica_stepDown(:,3);
    silica_stepDowni5f0p5 = silica_stepDown(:,4);
    silica_stepDowni5f0p1 = silica_stepDown(:,5);

    %% Solution
    i1f1 = stepShear(obj, i1, f1, silica_stepDownTime);
    i2f2 = stepShear(obj, i2, f2, silica_stepDownTime);
    i3f3 = stepShear(obj, i3, f3, silica_stepDownTime);
    i4f4 = stepShear(obj, i4, f4, silica_stepDownTime);

    transient_error = (sum((norm(i1f1.stress-silica_stepDowni5f2p5))./mean(silica_stepDowni5f2p5)) + ...
                       sum((norm(i2f2.stress-silica_stepDowni5f1p0))./mean(silica_stepDowni5f1p0)) + ...
                       sum((norm(i3f3.stress-silica_stepDowni5f0p5))./mean(silica_stepDowni5f0p5)) + ...
                       sum((norm(i4f4.stress-silica_stepDowni5f0p1))./mean(silica_stepDowni5f0p1)))  ...
                     ./length(silica_stepDownTime);

%% Error
fObj = SS_error + transient_error;
end