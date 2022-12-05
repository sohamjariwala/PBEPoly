function dX = elasticStrain(obj, t, gamma_e, shearRate, logintMu)
gamma_max = obj.gamma_e_max(logintMu);
    
% dX = shearRate*(1 - obj.structure_shear_rate(shearRate,logintMu)/shearRate*(gamma_e/gamma_max)^4);

phi = obj.phi_a(logintMu);
d_f = obj.par.d_f;

Mu = exp(logintMu);
mu0 = Mu(1);

dMu = obj.momicDerivative5(t, logintMu, shearRate);
dMu0dt = dMu(1)*mu0;

dgamma_maxdt = -obj.cnst.phi_p*obj.gamma_lin*2/d_f*((phi - obj.phi_pc)).^((d_f-1)/(3-d_f))...
    ./(obj.phi_max(logintMu)- obj.phi_pc).^((2)/(3-d_f))*obj.cnst.a_p*mu0^(-3/d_f) ...
   .*dMu0dt;

gamma_dot_p = obj.gamma_dot_p(gamma_e, shearRate, logintMu);

dX = (gamma_dot_p - gamma_e/gamma_max*abs(gamma_dot_p))*(dgamma_maxdt >= 0) ...
   + (gamma_dot_p - gamma_e/gamma_max*abs(gamma_dot_p) + gamma_e/gamma_max*dgamma_maxdt)*(dgamma_maxdt<0);

% dX = (gamma_dot_p - gamma_e/gamma_max*abs(gamma_dot_p));

end