function out = stepShear(obj, initialShearRate, finalShearRate, time, init)
    % number of moments to be solved for
    numMoments = 5;

    % Step-up and Step-down
    shear_rate = @(t) finalShearRate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Obtaining initial conditions
    if nargin < 5
    % Value at initial steady state
    init = obj.steadyShearODE(initialShearRate);
    init.gamma_e = obj.gamma_e_max(init.logintMu);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Solving the system of ODEs
    % Setting tolerance for ODEs
    tstart = tic;
    odeopts = odeset('RelTol',1e-4,'AbsTol', 1e-4, 'Stats','off', 'Events',@(t,X) myEvent(t,X,tstart));
    % Solving for transient state at specified times
    fun = @(t, X) [obj.momicDerivative5(t, X(1:numMoments)', obj.gamma_dot_p(X(end-1), shear_rate(t), X(1:numMoments)')); 
                   obj.elasticStrain(t, X(end-1), shear_rate(t), X(1:numMoments)');
                   obj.shearStressDE(t, X(end), X(end-1), shear_rate(t), X(1:numMoments)')];
try
              
    out.sol = ode15s(fun, ...
        [0, time(end)], ...
        [init.logintMu, init.gamma_e, init.stress],...
        odeopts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Running the solution

    % Obtaining solution variables
    if length(time) > 2
        [~, msgid] = lastwarn;
            if strcmp(msgid, 'MATLAB:illConditionedMatrix') || strcmp(msgid, 'MATLAB:ode15s:IntegrationTolNotMet')
                error('Matrix is singular, close to singular or badly scaled. Results may be inaccurate')
            end
        lastwarn('')
        
        [sol_Eval, derivatives] = deval(out.sol, time);

        out.logintMu = sol_Eval(1:numMoments,:)'; 
        out.gamma_e = sol_Eval(end-1,:)';
        out.gamma_e_dot = derivatives(end-1,:)';
        for i = 1:length(time)
            out.phi_a(i) = obj.phi_a(out.logintMu(i,:));
            out.stress(i) = sol_Eval(end,i)';
        end

    else
        [~, msgid] = lastwarn;
            if strcmp(msgid, 'MATLAB:illConditionedMatrix') || strcmp(msgid, 'MATLAB:ode15s:IntegrationTolNotMet')
                error('Matrix is singular, close to singular or badly scaled. Results may be inaccurate')
            end
        lastwarn('')
        
        [sol_Eval, derivatives] = deval(out.sol, out.sol.x);
        time = out.sol.x;
        out.time = time;
        out.logintMu = sol_Eval(1:numMoments,:)'; 
        out.gamma_e = sol_Eval(end-1,:)';
        out.gamma_e_dot = derivatives(end-1,:)';
        for i = 1:length(time)
            out.phi_a(i) = obj.phi_a(out.logintMu(i,:));
            out.stress(i) = sol_Eval(end,i)';
        end      
    end

catch
    warning('Problem evaluating ODE:STEPSHEAR');
    out.phi_a = 10^6*ones(size(time));
    out.logintMu = init.logintMu*ones(numMoments,length(time));
    out.gamma_e = 10^6*ones(size(time));
    lastwarn('')
    out.EXITFLAG = -10;
    out.stress = 10^6*ones(size(time));
    out.elastic_comp = 10^6*ones(size(time));
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output variables
    for i = 1:length(time)
    % Computing overall macroscopic stress
    out.elastic_comp(i) = obj.elastic_stress(out.logintMu(i,:),out.gamma_e(i));
    end
    out.elastic_comp = out.elastic_comp';
end

