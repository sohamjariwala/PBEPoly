function out = UDLAOS(obj, gamma_0, omega, time, init)
% UD-LAOS stress calculation
numMoments = 5;
% Sinusoidal shear rate
shear_rate = @(t) gamma_0*omega + gamma_0*omega*cos(omega*t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Obtaining initial conditions
    if nargin < 5
    % Value at initial steady state
    init = obj.steadyShear(gamma_0*omega);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Solving the system of ODEs
    % Setting tolerance for ODEs
    tstart = tic;
    odeopts = odeset('RelTol',1e-6,'AbsTol',1e-6,'Stats','on','Events',@(t,X) obj.myEvent(t,X,tstart));
    % Solving for transient state at specified times
    fun = @(t, X) [obj.momicDerivative5(t, X(1:numMoments)', shear_rate(t)); 
                   obj.elasticStrain(t, X(end), shear_rate(t), X(1:numMoments)')];
% try
    out.sol = ode15s(fun, ...
        [0, time(end)], ...
        [init.logintMu, init.gamma_e],...%
        odeopts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Running the solution

    % Obtaining solution variables
    if length(time) > 2
        [sol_Eval, derivatives] = deval(out.sol, time);

        out.logintMu = sol_Eval(1:numMoments,:)'; 
        out.phi_a = obj.phi_a(out.logintMu);
        out.gamma_e = sol_Eval(end,:)';
        out.gamma_e_dot = derivatives(end,:)';

        [~, msgid] = lastwarn;
            if strcmp(msgid, 'MATLAB:illConditionedMatrix') || strcmp(msgid, 'MATLAB:ode15s:IntegrationTolNotMet')
                error('Matrix is singular, close to singular or badly scaled. Results may be inaccurate')
            end
        lastwarn('')
    else
        [sol_Eval, derivatives] = deval(out.sol, out.sol.x);
        time = out.sol.x;
        out.time = time;
        out.logintMu = sol_Eval(1:numMoments,:)'; 
        out.phi_a = obj.phi_a(out.logintMu);
        out.gamma_e = sol_Eval(end,:)';
        out.gamma_e_dot = derivatives(end,:)';
        
    end
    
% % catch
%     warning('Problem evaluating ODE:UDLAOS');
%     out.phi_a = 10^6*ones(size(time));
%     out.logintMu = init.logintMu*ones(numMoments,length(time));
%     out.gamma_e = 10^6*ones(size(time));
%     lastwarn('')
%     out.EXITFLAG = -10;
%     out.stress = 10^6*ones(size(time));
%     out.elastic_comp = 10^6*ones(size(time));
%     return
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output variables
    for i = 1:length(time)
    % Computing overall macroscopic stress
    out.stress(i) = obj.total_stress(out.logintMu(i,:), ...
                                  out.gamma_e(i), ...
                                  shear_rate(time(i)));
    out.elastic_comp(i) = obj.elastic_stress(out.logintMu(i,:),out.gamma_e(i));
    end
    out.stress = out.stress';
    out.elastic_comp = out.elastic_comp';
end