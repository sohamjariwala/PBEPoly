addpath('../');
format long
rng default

%% Starting point
initialPoint = [   0.010453871382087   0.849999156000000  -5.060780142936760   2.119061369000000   0.957866855000000   6.254145353393902 0];
exploreRange = 0.10;
lb = min([(1-exploreRange)*initialPoint; (1+exploreRange)*initialPoint]);
ub = max([(1-exploreRange)*initialPoint; (1+exploreRange)*initialPoint]);

lb_max = [ 0.008, 0.40, -9, 2.01, 0.5, 4.6, 0.1];
ub_max = [5, 0.99, -3, 2.6, 0.99, 6.6, 11];
lb = max([lb;lb_max]);
ub = min([ub;ub_max]);

% Initializing the fluid object and functions with default parameter
fluid = PBEPoly;

%% Create the file for parameters

writematrix(["N", "W", "alfa", "b_0", "d_f", "porosity", "m_p", ...
    "kh", "error"], ...
    "Fits_ArmstrongSilica.csv", ...
    "WriteMode","overwrite", ...
    "Delimiter",",");

%% Fit loop
ITER=100; Nruns = 15;
parfor run = 1:Nruns
    ERBEST=100;
    parVec = initialPoint;
    parVecBest = parVec;
    obj = fluid;

    for ijk = 1:ITER
    
        fprintf("Iteration number #%d\n",ijk );

        if rem(ijk,10)==0
            fprintf("Current BEST PARAMETER values = \n")
            fprintf("%1.8f ",parVec);
            fprintf("\n");
        end

        if ijk<=1/2*ITER
            MULT=.15;
        elseif ijk>1/2*ITER && ijk <=3/4*ITER
            MULT=.075;
        else
            MULT=.025;
        end

        fprintf("Current PARAMETER values = \n")
        parVecNew = sign(parVec).*(sqrt(abs(parVec) + MULT*(0.5-rand(size(parVec))).*((ub-lb)))).^2;
        fprintf("%1.8f ",parVecNew);
        fprintf("\n");

        % Objective function
        func = @(x) objectiveFunctionSilica(obj, x);

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

    % Assigning values to the parameters
    obj.par.W = exp(parVecBest(1))-1;
    obj.par.alfa = parVecBest(2);
    obj.par.b_0 = exp(parVecBest(3));
    obj.par.d_f = parVecBest(4);
    obj.par.porosity = parVecBest(5);
    obj.par.m_p = exp(parVecBest(6));
    obj.par.kh = parVecBest(7);
    
    disp(obj.par)

    N = run;
    W = obj.par.W;
    alfa =  obj.par.alfa;
    b_0 =  obj.par.b_0;
    d_f =  obj.par.d_f;
    porosity =  obj.par.porosity;
    m_p =  obj.par.m_p;
    kh = obj.par.kh;
    error = ERBEST;

    % Write to file
    writematrix([N, W, alfa, b_0, d_f, porosity, m_p, ...
    kh, error], ...
    "Fits_ArmstrongSilica.csv", ...
    "Delimiter",",", ...
    "WriteMode","append");
end