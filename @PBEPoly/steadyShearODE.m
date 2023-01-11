function out = steadyShearODE(obj, shear_rate, initialConditions)
% Steady shear stress calculation

if nargin > 2
    init.logintMu=initialConditions.logintMu;
    init.stress = initialConditions.stress;
    init.gamma_e = initialConditions.gamma_e;
else
    init.logintMu=obj.InitialCondition;
    init.gamma_e = obj.gamma_lin;
    init.stress = 100;
end

numMoments = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Solving the system of ODEs
    % Setting tolerance for ODEs
    tstart = tic;
    odeopts = odeset('RelTol',1e-4,'AbsTol', 1e-4, 'Stats','off', 'Events',@(t,X) myEvent(t,X,tstart));
    % Solving for transient state at specified times
    fun = @(t, X) [obj.momicDerivative5(t, X(1:numMoments)', obj.gamma_dot_p(X(end-1), shear_rate, X(1:numMoments)')); 
                   obj.elasticStrain(t, X(end-1), shear_rate, X(1:numMoments)');
                   obj.shearStressDE(t, X(end), X(end-1), shear_rate, X(1:numMoments)')];
try
              
    out.sol = ode15s(fun, ...
        [0, 1e5], ...
        [init.logintMu, init.gamma_e, init.stress],...
        odeopts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Running the solution

        [sol_Eval, derivatives] = deval(out.sol, out.sol.x(end));

        out.logintMu = sol_Eval(1:numMoments,end)'; 
        out.gamma_e = sol_Eval(end-1,end)';
        out.gamma_e_dot = derivatives(end-1,end)';
        out.phi_a = obj.phi_a(out.logintMu(end,:));
        out.stress = sol_Eval(end,end)';
        [~, msgid] = lastwarn;
            if strcmp(msgid, 'MATLAB:illConditionedMatrix') || strcmp(msgid, 'MATLAB:ode15s:IntegrationTolNotMet')
                error('Matrix is singular, close to singular or badly scaled. Results may be inaccurate')
            end
        lastwarn('')
        out.EXITFLAG = 1;
    
catch
    warning('Problem evaluating ODE:STEPSHEAR');
    out.phi_a = 10^6*ones(1);
    out.logintMu = 10^6*init.logintMu*ones(1);
    out.gamma_e = 10^6*ones(1);
    lastwarn('')
    out.EXITFLAG = -10;
    out.stress = 10^6*ones(1);
    out.elastic_comp = 10^6*ones(1);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output variables

    out.elastic_comp = obj.elastic_stress(out.logintMu,out.gamma_e);
end

