addpath('../');
format long
rng default
%% Starting point
par = [0.00533651120435152,0.97,3.87795454688471e-05,2.37042495495754,...
       1,324.497448647856,1085.52243880163,...
       10.582149330422,0.0461599694066823,0.0227165642151505];

initialPoint = [log(par(1)+1), par(2), log(par(3)), par(4), par(5), log(par(6)), par(7), par(8), par(9), par(10)];

exploreRange = 0.05;
lb = min([(1-exploreRange)*initialPoint; (1+exploreRange)*initialPoint]);
ub = max([(1-exploreRange)*initialPoint; (1+exploreRange)*initialPoint]);

% Initializing the fluid object and functions with default parameter
fluid = PBEPoly;

%% Create the file for parameters

writematrix(["N", "W", "alfa", "b_0", "d_f", ...
    "porosity", "m_p", "G_0", "sigma_y0", "phi_p" ...
    "mu_s", "error"], ...
    "Fits_Hippetal.csv", ...
    "WriteMode","overwrite", ...
    "Delimiter",",");

%% Fit loop
ITER=100; Nruns = 15;
parfor run = 1:Nruns
    ERBEST=100;
    parVec = initialPoint;
    parVecBest = parVec;
    obj = fluid;
    obj.par.kh = 0;
    obj.cnst.a_p = 10e-9;
    obj.par.p = 0;


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
        func = @(x) objectiveFunctionCarbonBlack(obj, x);

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
    obj.cnst.G_0 = parVecBest(7);
    obj.cnst.sigma_y0 = parVecBest(8);
    obj.cnst.phi_p = parVecBest(9);
    obj.cnst.mu_s = parVecBest(10);
    obj.par.kh = 0;


    disp(obj.par)

    N = run;
    W = obj.par.W;
    alfa =  obj.par.alfa;
    b_0 =  obj.par.b_0;
    d_f =  obj.par.d_f;
    porosity =  obj.par.porosity;
    m_p =  obj.par.m_p;
    G_0 = obj.cnst.G_0;
    sigma_y0 = obj.cnst.sigma_y0;
    phi_p = obj.cnst.phi_p;
    mu_s = obj.cnst.mu_s;
    error = ERBEST;

    % Write to file
    writematrix([N, W, alfa, b_0, d_f, porosity, m_p, G_0, sigma_y0, phi_p...
    mu_s, error], ...
    "Fits_Hippetal.csv", ...
    "Delimiter",",", ...
    "WriteMode","append");
end