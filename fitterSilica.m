%% Starting point
load('silicaTest4.mat')
    parVec(1) = log(1+silica.par.W);
    parVec(2) = silica.par.alfa;
    parVec(3) = log(silica.par.b_0*0.5);
    parVec(4) = silica.par.porosity;


lb = min([0.85*parVec; 1.15*parVec]);
ub = max([0.85*parVec; 1.15*parVec]);

ITER=1000;   ERBEST=100; Nruns = 15;
for run = 1:Nruns
    for ijk = 1:ITER
    
        fprintf("Run #%d | Iteration number #%d\n", run, ijk );

        if rem(ijk,10)==0
            fprintf("Current best parameter values = \n")
            disp(parVec);
            fprintf("Error = \n");
            disp(ERBEST);
        end

        if ijk<=1/2*ITER
            MULT=.15;
        elseif ijk>1/2*ITER & ijk <=3/4*ITER
            MULT=.075;
        else
            MULT=.025;
        end
%         fprintf("Current parameter values = \n")
        parVecNew = sign(parVec).*(sqrt(abs(parVec) + MULT*(0.5-rand(size(parVec))).*((ub-lb)))).^2;
%         disp(parVecNew);

        % Initializing the fluid object and functions with default parameter
        fluid(run) = PBEPoly;
        % Changing constants and parameters to ones obtained from monodisperse 
        % solution changed parameters
        fluid(run).cnst.phi_p = silica.cnst.phi_p;
        fluid(run).cnst.a_p = silica.cnst.a_p;
        fluid(run).cnst.mu_s = silica.cnst.mu_s;
        fluid(run).cnst.G_0 = silica.cnst.G_0;
        fluid(run).cnst.sigma_y0 = silica.cnst.sigma_y0;

        fluid(run).par.W = silica.par.W;
        fluid(run).par.alfa = silica.par.alfa;
        fluid(run).par.b_0 = silica.par.b_0;
        fluid(run).par.d_f = silica.par.d_f;
        fluid(run).par.porosity = silica.par.porosity;
        fluid(run).par.m_p = silica.par.m_p;

        obj = fluid(run);
        % Objective function
        func = @(x) objectiveFunctionSilica(obj, x);

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
    fluid(run).par.porosity = x(4);

    
    N = run;
    W = fluid(run).par.W;
    alfa =  fluid(run).par.alfa;
    b_0 =  fluid(run).par.b_0;
    d_f =  fluid(run).par.d_f;
    porosity =  fluid(run).par.porosity;
    m_p =  fluid(run).par.m_p;
    error = ERBEST;
    
    if run ==1
    % Writing values to file
    fprintf("Writing the values to file\n");
    TT = table(N, W, alfa, b_0, porosity, error);

    writetable(TT,...
    "Fits_Silica.csv",...
    'WriteVariableNames',1)
    else
       T2 = readtable('Fits_Silica.csv');
       T3 = table(N, W, alfa, b_0, porosity, error);
       TT = [T2;T3];
     writetable(TT,...
     "Fits_Silica.csv",...
    'WriteVariableNames',1)
    end

end