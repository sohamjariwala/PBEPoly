function out = steadyShear(obj, shear_rate, initialConditions)
% Steady shear stress calculation
 try
    % Initial guesses
    if nargin > 2
        init.logintMu= initialConditions.logintMu;
        init.stress =  initialConditions.stress;
    else
        init.logintMu=obj.InitialCondition;
        init.stress = 10;
    end

    % Solving for steady state using MATLAB fsolve
    fun = @(x) obj.momicDerivative5(0, x, shear_rate);
    options = optimoptions('fsolve','TolFun', eps, 'Display','off', 'FunctionTolerance', eps);
    [out.logintMu,~,out.EXITFLAG,~,~] = fsolve(fun, init.logintMu, options);
    
    % Output variables
    out.stress = total_stress_SS(obj, out.logintMu, shear_rate);
    out.sigma_y = obj.sigma_y(out.logintMu);

 catch
    warning('FSOLVE returned bad solution')
    out.stress = 10^6*ones(size(shear_rate));
    out.phi_a = 10^6*ones(size(shear_rate));
    out.gamma_e = 10^6*ones(size(shear_rate));
    out.logintMu=zeros(length(shear_rate), 5);
    out.EXITFLAG = -10;
 end
 
end