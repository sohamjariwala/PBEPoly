addpath('../');
format long
%% Starting point
initialPoint = [0.008, 0.62, -7.1062, 2.11, 0.897, 6.1485];
exploreRange = 0.20;

lb = min([(1-exploreRange)*initialPoint; (1+exploreRange)*initialPoint]);
ub = max([(1-exploreRange)*initialPoint; (1+exploreRange)*initialPoint]);

% Initializing the fluid object and functions with default parameter
fluid = PBEPoly;

%% Create the file for parameters

writematrix(["N", "W", "alfa", "b_0", "d_f", "porosity", "m_p", ...
    "G_0", "sigma_y0","mu_s", "error"], ...
    "Fits_ArmstrongSilica.csv", ...
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
            fprintf("Current BEST PARAMETER values = \n")
            fprintf("%1.8f ",parVec);
            fprintf("\n");
        end

        if ijk<=1/2*ITER
            MULT=.15;
        elseif ijk>1/2*ITER & ijk <=3/4*ITER
            MULT=.075;
        else
            MULT=.025;
        end

        fprintf("Current PARAMETER values = \n")
        parVecNew = sign(parVec).*(sqrt(abs(parVec) + MULT*(0.5-rand(size(parVec))).*((ub-lb)))).^2;
        fprintf("%1.8f ",parVecNew);
        fprintf("\n");

        obj = fluid;

        % Objective function
        func = @(x) objectiveFunctionSilica(fluid, x);

        ERROR_TOT = func(parVecNew);
        fprintf("Current ERROR value = \n")
        fprintf("%1.8f ",ERROR_TOT);
        fprintf("\n");
        
        if ERROR_TOT<=ERBEST
            parVecBest = parVecNew;
            parVec = parVecNew;
            ERBEST = ERROR_TOT;
        end
        fprintf("Current BEST ERROR value = \n")
        fprintf("%1.8f ",ERBEST);
        fprintf("\n");
    end
    x = parVecBest;

    % Assigning values to the parameters
    obj.par.W = exp(x(1))-1;
    obj.par.alfa = x(2);
    obj.par.b_0 = exp(x(3));
    obj.par.d_f = x(4);
    obj.par.porosity = x(5);
    obj.par.m_p = exp(x(6));
    
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
    "Fits_ArmstrongSilica.csv", ...
    "Delimiter",",", ...
    "WriteMode","append");
end