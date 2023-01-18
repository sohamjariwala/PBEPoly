function out = steadyShearODE(obj, shear_rate, initialConditions)
% Steady shear stress calculation using ODE solution at steady state

if nargin > 2
    init.logintMu = initialConditions.logintMu;
    init.stress = initialConditions.stress;
    init.A = initialConditions.A;
else
    init.logintMu=obj.InitialCondition;
    init.stress = 300;
    init.A = 1;
end

numMoments = 5;
    
%% Solving the system of ODEs
% Setting tolerance for ODEs
tstart = tic;
odeopts = odeset('RelTol',1e-2,'AbsTol',1e-4,'Stats','off','Events', ...
    @(t,X) obj.odeEvent(t,X,tstart));

% Simultaneous ODEs to be solved till stationary state is attained
fun = @(t, X) [obj.momicDerivative5(t, X(1:numMoments)', obj.gamma_dot_p(X(end), X(end-1), X(1:numMoments)',shear_rate));
    obj.Adot(t,X(end-1), X(1:numMoments)',X(end),shear_rate);
    obj.shearStressDE(t, X(end), shear_rate, X(1:numMoments)',X(end-1))];

try
    %% Running the solution
    out.sol = ode15s(fun, ...
        [0, 1e14*obj.tau(init.logintMu)], ...
        [init.logintMu, init.A, init.stress],...
        odeopts);

    [out.sol_Eval, out.derivatives] = deval(out.sol, out.sol.x(end));

    out.logintMu = out.sol_Eval(1:numMoments,end)';
    out.phi_a = obj.phi_a(out.logintMu(end,:));
    out.sigma_y = obj.sigma_y(out.logintMu(end,:));
    out.A = out.sol_Eval(end-1,end);
    out.stress = out.sol_Eval(end,end);
    out.viscous_stress = obj.viscous_stress(out.logintMu(end,:),shear_rate);
    out.back_stress = out.stress - out.viscous_stress - out.sigma_y; 
    out.EXITFLAG = 1;
    
catch
    warning('Problem evaluating ODE:STEADYSHEARODE');
    out.logintMu = 10^6*init.logintMu*ones(1);
    out.phi_a = 10^6*ones(1);
    out.sigma_y = 10^6*ones(1);
    out.A = 10^6*ones(1);
    out.stress = 10^6*ones(1);
    out.EXITFLAG = -10;
    lastwarn('')
    return
end
end