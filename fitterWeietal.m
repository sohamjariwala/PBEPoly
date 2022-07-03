%% Starting point
% parVec = [3.018240e+00, 5.848876e-01, -6.412818e+00, 2.272982e+00, 8.816055e-01, 5.944847e+00, 3.712128e+01, 7.804288e-01, 3.360543e-01];
parVec = [3.013704391154364   0.590198260682092  -6.354991898429272   2.268634081422561   0.883034591106096   5.953417691287351  37.277275082280944   0.782275208105718 0.333520407166337];

parVecNew = parVec;

lb = min([0.975*parVec; 1.025*parVec]);
ub = max([0.975*parVec; 1.025*parVec]);

ITER=100;   ERBEST=100; Nruns = 15;
for run = 1:Nruns
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

        % Initializing the fluid object and functions with default parameter
        fluid(run) = PBEPoly;
        % Changing constants and parameters to ones obtained from monodisperse 
        % solution changed parameters
        fluid(run).cnst.phi_p = 0.0300;
        fluid(run).cnst.a_p = 8e-9;
        fluid(run).cnst.mu_s = 0.3592;
        fluid(run).cnst.G_0 = 2.3896;
        fluid(run).cnst.sigma_y0 = 0.7951;

        fluid(run).par.W = 14.1442;
        fluid(run).par.alfa = 0.4096;
        fluid(run).par.b_0 = 0.0015;
        fluid(run).par.d_f = 2.3325;
        fluid(run).par.porosity = 0.8192;
        fluid(run).par.m_p = 978;

        obj = fluid(run);
        % Objective function
        func = @(x) objectiveFunctionWei(obj, x);

        ERROR_TOT = func(parVecNew);

        if ERROR_TOT<ERBEST
            parVecBest = parVecNew;
            parVec = parVecNew;
            ERBEST = ERROR_TOT;
        end

    end
    x = parVecBest;
    % Assigning values to the parameters
    fluid(run).par.W = exp(x(1))-1;
    fluid(run).par.alfa = x(2);
    fluid(run).par.b_0 = exp(x(3));
    fluid(run).par.d_f = x(4);
    fluid(run).par.porosity = x(5);
    fluid(run).par.m_p = exp(x(6));
    fluid(run).cnst.G_0 = x(7);
    fluid(run).cnst.sigma_y0 = x(8);
    fluid(run).cnst.mu_s = x(9);
    
    N = run;
    W = fluid(run).par.W;
    alfa =  fluid(run).par.alfa;
    b_0 =  fluid(run).par.b_0;
    d_f =  fluid(run).par.d_f;
    porosity =  fluid(run).par.porosity;
    m_p =  fluid(run).par.m_p;
    G_0 = fluid(run).cnst.G_0;
    sigma_y0 = fluid(run).cnst.sigma_y0;
    mu_s = fluid(run).cnst.mu_s;
    error = ERBEST;
    
    if run ==1
    % Writing values to file
    fprintf("Writing the values to file\n");
    TT = table(N, W, alfa, b_0, d_f, porosity, m_p, G_0, sigma_y0,...
                   mu_s, error);

    writetable(TT,...
    "Fits_Weietal.csv",...
    'WriteVariableNames',1)
    else
       T2 = readtable('Fits_Weietal.csv');
       T3 = table(N, W, alfa, b_0, d_f, porosity, m_p, G_0, sigma_y0,...
                   mu_s, error);
       TT = [T2;T3];
     writetable(TT,...
     "Fits_Weietal.csv",...
    'WriteVariableNames',1)
    end

end