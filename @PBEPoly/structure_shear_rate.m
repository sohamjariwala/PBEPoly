function x = structure_shear_rate(obj,shearRate,logintMu)
phi_a = obj.phi_a(logintMu);

% Shear rate corresponding to steady shear structure 
A =  - obj.par.b_0 * obj.cnst.phi_p ...
    * ((phi_a/obj.cnst.phi_p)^((2)/(obj.par.d_f - 3)) ...
    -  obj.par.m_p^(1/obj.par.d_f)*(phi_a/obj.cnst.phi_p)^((3)/(obj.par.d_f - 3)));

B = 4 * cutOff(obj, logintMu) * obj.par.alfa ...
    * (obj.cnst.phi_p^2/pi) ...
    * (phi_a/obj.cnst.phi_p)^((obj.par.d_f)/(obj.par.d_f - 3));

% eta = etaTrimodal(obj, logintMu);
phi_a = obj.phi_a(logintMu);
eta = obj.etaTrimodal(logintMu);

C = 2 * cutOff(obj, logintMu) * obj.cnst.k_b * obj.cnst.T * obj.cnst.phi_p^2 ...
   ./( 2 * obj.cnst.mu_s * eta * obj.par.W * pi * obj.cnst.a_p^3) ...
   * (phi_a/obj.cnst.phi_p).^((obj.par.d_f + 3)/(obj.par.d_f - 3));

structure_SR = abs((B + sqrt(B.^2-4.*A.*C))./(2.*A));
% x = structure_SR;
x = max(shearRate, min(1e2,structure_SR));
% structure_SR <= shearRate
end