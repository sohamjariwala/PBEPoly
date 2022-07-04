function dX = shearStressDE(obj, t, sigma, gamma_e, shearRate, logintMu)
% Differential equation for shear stress based on Saramito model with
% explicit yield stress and thermodynamically consistent
% elasto-visco-plasticity.
phi = obj.phi_a(logintMu);

dX = - obj.G(logintMu) * ...
    max(0, abs(sigma - obj.cnst.sigma_y0/obj.gamma_lin*gamma_e)...
    /obj.eta(phi,logintMu))*sign(sigma) ...
    + obj.G(logintMu)*shearRate;

end