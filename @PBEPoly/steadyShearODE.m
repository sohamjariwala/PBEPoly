function out = steadyShearODE(obj, shear_rate, initialConditions)
% Steady shear stress calculation

if nargin > 2
    loginitialMu=initialConditions;
else
    loginitialMu=obj.InitialCondition;
end

% Setting tolerance for ODEs
odeopts = odeset('Event',@eventfun,'RelTol',1e-4,'AbsTol',1e-4,'Stats','off');

% Solving for transient state at specified times
fun = @(t, X) [obj.momicDerivative6(t, X(1:5)', shear_rate); 
               obj.elasticStrain(t, X(end), shear_rate, X(1:5)')];
           
out.sol = ode15s(fun, ...
            [0, 20], ...
            [loginitialMu, obj.gamma_lin],...
            odeopts);
 
        
function [x,isterm,dir] = eventfun(t,X)
    dy = fun(t, X);     
    x = norm(dy) - 1e-6;
    isterm = 1;
    dir = 0;  %or -1, doesn't matter
end
        
        
try
    % Obtaining solution variables
    [sol_Eval, derivatives] = deval(out.sol, out.sol.x(end));
    
    out.logintMu = sol_Eval(1:6,end)'; 
    out.phi_a = obj.phi_a(out.logintMu);
    out.gamma_e =  sol_Eval(end,end)';
    out.gamma_e_dot = derivatives(end,end)';

    [~, msgid] = lastwarn;
    if strcmp(msgid, 'MATLAB:illConditionedMatrix')
        out.EXITFLAG = -10;
        error('Matrix is singular, close to singular or badly scaled. Results may be inaccurate');
    else
        out.EXITFLAG = 1;
    end
    lastwarn('')
    % Computing overall macroscopic stress
    out.stress = obj.total_stress(out.logintMu, ...
                                  obj.gamma_lin, ...
                                  shear_rate);
catch
    warning('Problem evaluating ODE:UDLAOS');
    out.phi_a = 10^6*ones(1);
    out.gamma_e = 10^6*ones(1);
    lastwarn('')
    out.EXITFLAG = -10;
    out.stress = 10^6;
    out.logintMu = obj.InitialCondition;
end
end

