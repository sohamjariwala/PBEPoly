addpath('../');
%% Starting point
initialPoint = [3.125965752659946   0.593094235000000  -5.770599429453542   2.273852404000000   0.885542089000000   5.686169744124814  36.288052170000000   0.819958230000000   0.325495437000000];
parVecNew = initialPoint;
exploreRange = 0.10;

lb = min([(1-exploreRange)*initialPoint; (1+exploreRange)*initialPoint]);
ub = max([(1-exploreRange)*initialPoint; (1+exploreRange)*initialPoint]);

% Initializing the fluid object and functions with default parameter
fluid = PBEPoly;
% Changing constants and parameters to ones obtained from monodisperse 
% solution changed parameters
fluid.cnst.phi_p = 0.0300;
fluid.cnst.a_p = 8e-9;
fluid.cnst.mu_s = 0.3592;
fluid.cnst.G_0 = 2.3896;
fluid.cnst.sigma_y0 = 0.7951;

fluid.par.W = 14.1442;
fluid.par.alfa = 0.4096;
fluid.par.b_0 = 0.0015;
fluid.par.d_f = 2.3325;
fluid.par.porosity = 0.8192;
fluid.par.m_p = 978;
fluid.par.p = 4;


%% Create the file for parameters

writematrix(["N", "W", "alfa", "b_0", "d_f", "porosity", "m_p", ...
    "G_0", "sigma_y0","mu_s", "error"], ...
    "Fits_Weietal.csv", ...
    "WriteMode","overwrite", ...
    "Delimiter",",");

%% Fit loop
ITER=100; Nruns = 15;
parfor run = 1:Nruns
    ERBEST=100;
    parVec = initialPoint;
    parVecBest = parVec;

    for ijk = 1:ITER
    
        fprintf("Iteration number #%d\n",ijk );

        if rem(ijk,10)==0
            fprintf("Current best parameter values = \n")
            disp(parVec);
        end

        if ijk<=1/2*ITER
            MULT=.15;
        elseif ijk>1/2*ITER & ijk <=3/4*ITER
            MULT=.075;
        else
            MULT=.025;
        end
        fprintf("Current parameter values = \n")
        parVecNew = sign(parVec).*(sqrt(abs(parVec) + MULT*(0.5-rand(size(parVec))).*((ub-lb)))).^2;
        disp(parVecNew);

        obj = fluid;
        % Objective function
        func = @(x) objectiveFunctionWei(obj, x);

        ERROR_TOT = func(parVecNew);

        if ERROR_TOT<=ERBEST
            parVecBest = parVecNew;
            parVec = parVecNew;
            ERBEST = ERROR_TOT;
        end

    end
    x = parVecBest;

    % Assigning values to the parameters
    obj.par.W = exp(x(1))-1;
    obj.par.alfa = x(2);
    obj.par.b_0 = exp(x(3));
    obj.par.d_f = x(4);
    obj.par.porosity = x(5);
    obj.par.m_p = exp(x(6));
    obj.cnst.G_0 = x(7);
    obj.cnst.sigma_y0 = x(8);
    obj.cnst.mu_s = x(9);
    
    N = run;
    W = obj.par.W;
    alfa =  obj.par.alfa;
    b_0 =  obj.par.b_0;
    d_f =  obj.par.d_f;
    porosity =  obj.par.porosity;
    m_p =  obj.par.m_p;
    G_0 = obj.cnst.G_0;
    sigma_y0 = obj.cnst.sigma_y0;
    mu_s = obj.cnst.mu_s;
    error = ERBEST;

    % Write to file
    writematrix([N, W, alfa, b_0, d_f, porosity, m_p, ...
    G_0, sigma_y0,mu_s, error], ...
    "Fits_Weietal.csv", ...
    "Delimiter",",", ...
    "WriteMode","append");
end