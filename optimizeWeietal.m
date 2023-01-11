% Script to optimize the values for parameters in the polydisperse
% population balance consitutive model.
clear; clc
Nfits = 15;

% Initial vector of parameters
numParameters = 9;

% Initial guess
% parVec = [2.975301805672655   0.600266541744971  -6.212522743705272   2.079956662255416   0.871154549640630   6.067781189331943  39.578240245118025   0.489379810247538   0.201637100515005];
% parVec = [2.905571634707300   0.585389562203110  -6.057648176214162   2.131929000429122   0.891764092345692   5.917187722441486  39.067565659463568   0.7951   0.196624904006037];
% parVec = [2.864550050745758   0.572115011942363  -5.941326101945536   2.33153231558815501 0.902963313001152   5.966540711533419  39.463950524306519   0.776901649455502   0.3592200855597641966];
% parVec = [2.936134665904657   0.560929744046274  -5.982037049863742   2.286001277752308   0.896793889319817   6.115461098553929  39.024089138217356   0.771603183329606   0.356958563860062];
% parVec = [3.009184927963146   0.574820100362079  -6.129916683018121   2.284651220139327   0.896239432060398   5.967239794876706  38.363935423635844   0.771529572427426   0.348133957695984];
% parVec = [3.024952387572246   0.589041528128687  -6.249040808649854   2.270036781987128   0.887832090476037   5.840684821201803  38.501763148417659   0.790236543136271   0.339578696109709];
% parVec = [3.058294210350837   0.590840162774347  -6.279278825175675   2.265820000365196   0.883692869941822   5.888818412949577  38.135980833025492   0.790651168745462   0.336344273576225];
% parVec = [3.052969e+00, 5.902315e-01, -6.315614e+00, 2.265117e+00, 8.821384e-01, 5.878319e+00, 3.806119e+01, 7.888523e-01, 3.362045e-01];
parVec = [3.018240e+00, 5.848876e-01, -6.412818e+00, 2.272982e+00, 8.816055e-01, 5.944847e+00, 3.712128e+01, 7.804288e-01, 3.360543e-01];

lb = min([0.975*parVec; 1.025*parVec]);
ub = max([0.975*parVec; 1.025*parVec]);

% Fitting using Genetic algorithm
W = [];
alfa = W; b_0 = W; d_f = W;
porosity = W; m_p = W;
error = W; N = W; G_0 = W; sigma_y0 = W; mu_s = W;

% i =1;
for i = 1:Nfits
fprintf("Starting Run %d\n", i);

% Initializing the fluid object and functions with default parameter
fluid(i) = PBEPoly;
% Changing constants and parameters to ones obtained from monodisperse 
% solution changed parameters
fluid(i).cnst.phi_p = 0.0300;
fluid(i).cnst.a_p = 8e-9;
fluid(i).cnst.mu_s = 0.3592;
fluid(i).cnst.G_0 = 2.3896;
fluid(i).cnst.sigma_y0 = 0.7951;

fluid(i).par.W = 14.1442;
fluid(i).par.alfa = 0.4096;
fluid(i).par.b_0 = 0.0015;
fluid(i).par.d_f = 2.3325;
fluid(i).par.porosity = 0.8192;
fluid(i).par.m_p = 978;

obj = fluid(i);
% Objective function
func = @(x) objectiveFunctionWei(obj, x);


% Genetic algorithm
%     options = optimoptions('ga','UseParallel',true, 'Display', 'iter');
%     [x,fval,exitflag] = ga(func,numParameters,[],[],[],[],lb,ub,[],options);

% % Parallel tempering
%    x = ParTemp(func,lb+rand*(ub-lb),1e-2,lb,'MIN',ub,'MAX');

% % Simulated annealing
options = optimoptions('simulannealbnd', 'PlotFcns',...
          {@saplotbestx,@saplotbestf,@saplotx,@saplotf},...
          'OutputFcn', @myoutSA);
   [x,fval,exitflag,output]  = simulannealbnd(func, lb+rand*(ub-lb), lb, ub, options);
    
% Assigning values to the parameters
    fluid(i).par.W = exp(x(1))-1;
    fluid(i).par.alfa = x(2);
    fluid(i).par.b_0 = exp(x(3));
    fluid(i).par.d_f = x(4);
    fluid(i).par.porosity = x(5);
    fluid(i).par.m_p = exp(x(6));
    fluid(i).cnst.G_0 = x(7);
    fluid(i).cnst.sigma_y0 = x(8);
    fluid(i).cnst.mu_s = x(9);
    
    N = i;
    W = fluid(i).par.W;
    alfa =  fluid(i).par.alfa;
    b_0 =  fluid(i).par.b_0;
    d_f =  fluid(i).par.d_f;
    porosity =  fluid(i).par.porosity;
    m_p =  fluid(i).par.m_p;
    G_0 = fluid(i).cnst.G_0;
    sigma_y0 = fluid(i).cnst.sigma_y0;
    mu_s = fluid(i).cnst.mu_s;
    error = fval;
    
    if i ==1
    % Writing values to file
    fprintf("Writing the values to file\n");
    TT = table(N, W, alfa, b_0, d_f, porosity, m_p, G_0, sigma_y0,...
                   mu_s, error);

    writetable(TT,...
    "gaOptimized.csv",...
    'WriteVariableNames',1)
    else
       T2 = readtable('gaOptimized_Wei');
       T3 = table(N, W, alfa, b_0, d_f, porosity, m_p, G_0, sigma_y0,...
                   mu_s, error);
       TT = [T2;T3];
     writetable(TT,...
     "gaOptimized.csv",...
    'WriteVariableNames',1)
    end

end
