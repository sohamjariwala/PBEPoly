function dX = shearStressDE(obj, ~, sigma, shear_rate, logintMu, A)
% Differential equation for shear stress based on Saramito model with
% explicit yield stress and thermodynamically consistent
% elasto-visco-plasticity.
phi = obj.phi_a(logintMu);

dX = (...
    - sigma_eff(obj,sigma,logintMu, A) ...
    + obj.sigma_y(logintMu)/obj.structure_shear_rate(logintMu)*shear_rate ...
    + obj.cnst.mu_s*obj.eta(phi, logintMu)*shear_rate ...
    )...
    /obj.tau(logintMu);

end