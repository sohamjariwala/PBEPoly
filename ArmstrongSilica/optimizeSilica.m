addpath('../');
format long
rng default

%% Starting point
initialPoint = [0.08, 0.62, -7.1062, 2.11, 0.897, 6.1485, 0.7];

exploreRange = 0.50;
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
    "optim_ArmstrongSilica.csv", ...
    "WriteMode","overwrite", ...
    "Delimiter",",");

%% Fit loop
    obj = fluid;
    func = @(x) objectiveFunctionSilica(obj, x);
%     numParameters = length(lb);
%     options = optimoptions('ga','UseParallel',true, 'Display', 'iter', 'PlotFcn','gaplotbestf');
%     [x,fval,exitflag] = ga(func,numParameters,[],[],[],[],lb_max,ub_max,[],options);

%     options = optimoptions('simulannealbnd','PlotFcns',...
%           {@saplotbestx,@saplotbestf,@saplotx,@saplotf}, 'Display', 'iter');
%     x = simulannealbnd(func,initialPoint,lb_max,ub_max,options);


x = ParTemp(func,initialPoint,1,lb_max,'MIN',ub_max,'MAX');

    parVecBest = x;
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
    "optim_ArmstrongSilica.csv", ...
    "Delimiter",",", ...
    "WriteMode","append");