function dX = shearStressDE(obj, t, sigma, gamma_e, logintMu, gamma_dot)
% Differential equation for shear stress based on Saramito model with
% explicit yield stress and thermodynamically consistent
% elasto-visco-plasticity.

dX = - obj.G(logintMu) * ...
    max(0, abs(sigma - obj.par.sigma_y/obj.gamma_lin*gamma_e)...
    /obj.eta(logintMu))*sign(sigma) ...
    + obj.G(logintMu)*gamma_dot;

end