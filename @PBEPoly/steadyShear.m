function out = steadyShear(obj, shear_rate, initialConditions)
% Steady shear stress calculation
 try
    % Initial guesses
    if nargin > 2
        init.logintMu= initialConditions.logintMu;
    else
        init.logintMu=obj.InitialCondition;
    end

    % Solving for steady state using MATLAB fsolve
    fun = @(x) obj.momicDerivative5(0, x, obj.gamma_dot_p(obj.total_stress_SS(x, shear_rate), 1, x,shear_rate));
    options = optimoptions('fsolve','TolFun', eps, 'Display','off', 'FunctionTolerance', eps);
    [out.logintMu,out.derivatives,out.EXITFLAG,~,~] = fsolve(fun, init.logintMu, options);
    
    % Output variables
    out.stress = total_stress_SS(obj, out.logintMu, shear_rate);
    out.phi_a = obj.phi_a(out.logintMu(end,:));
    out.sigma_y = obj.sigma_y(out.logintMu);
    out.A = 1;

 catch
    warning('FSOLVE returned bad solution')
    out.stress = 10^6*ones(size(shear_rate));
    out.phi_a = 10^6*ones(size(shear_rate));
    out.logintMu = zeros(length(shear_rate), 5);
    out.sigma_y =  10^6*ones(size(shear_rate));
    out.EXITFLAG = -10*ones(size(shear_rate));
 end
 
end