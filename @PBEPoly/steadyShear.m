function out = steadyShear(obj, shear_rate)
% Steady shear stress calculation
 try
    % Solving for steady state using MATLAB fsolve
    fun = @(x) mwasameModelSS(obj, x, shear_rate);
    options = optimset('Display','off');

    [out.phi_a,~,out.EXITFLAG,~,~] = fsolve(fun, obj.phi_pc, options);

    % Output variables
    out.gamma_e = sqrt(shear_rate/...
        obj.structure_shear_rate(out.phi_a, shear_rate))*obj.gamma_lin;
    
    out.stress = total_stress(obj, out.phi_a, out.gamma_e, shear_rate);
    out.G = obj.G(out.phi_a);
 catch
    out.stress = zeros(size(shear_rate));
    out.phi_a = zeros(size(shear_rate));
    out.gamma_e = zeros(size(shear_rate));
    out.EXITFLAG = -10;
 end
end