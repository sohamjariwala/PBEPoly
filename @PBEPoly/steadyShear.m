function out = steadyShear(obj, shear_rate, initialConditions)
% Steady shear stress calculation
 try
    % Solving for steady state using MATLAB fsolve
    fun = @(x) momicDerivative5(obj, 0, x, shear_rate);
    options = optimoptions('fsolve','TolFun', eps, 'Display','off', 'FunctionTolerance', eps);
    
    if nargin > 2
        [out.logintMu,~,out.EXITFLAG,~,~] = fsolve(fun, initialConditions, options);
    else
        [out.logintMu,~,out.EXITFLAG,~,~] = fsolve(fun, obj.InitialCondition, options);
    end
    
    % Output variables
    out.gamma_e = obj.gamma_lin;
    out.stress = total_stress(obj, out.logintMu, out.gamma_e, shear_rate);
    out.G = obj.G(out.logintMu);
 catch
    out.stress = 10^6*ones(size(shear_rate));
    out.phi_a = 10^6*ones(size(shear_rate));
    out.gamma_e = 10^6*ones(size(shear_rate));
    out.logintMu=zeros(length(shear_rate), 5);
    out.EXITFLAG = -10;
 end
end