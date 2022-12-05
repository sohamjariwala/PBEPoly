function dX = shearStressDE(obj, t, sigma, gamma_e, shearRate, logintMu)
% Differential equation for shear stress based on Saramito model with
% explicit yield stress and thermodynamically consistent
% elasto-visco-plasticity.
phi = obj.phi_a(logintMu);

% dX = - obj.G(logintMu)*...
%     max(0, abs(sigma - obj.G(logintMu)/obj.gamma_lin*gamma_e^2)) ...
%     ./obj.eta(phi, logintMu)...
%     *sign(sigma) ...
%     + obj.G(logintMu)*shearRate;

dX = - max(0, abs(sigma - obj.G(logintMu)/obj.gamma_lin*gamma_e^2)) ...
    / obj.tau(logintMu,shearRate,gamma_e) * sign(sigma)...
    + (obj.cnst.mu_s*obj.eta(phi, logintMu))/obj.tau(logintMu,shearRate,gamma_e)*shearRate;
end