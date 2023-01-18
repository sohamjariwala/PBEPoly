function x = structure_shear_rate(obj,logintMu)
% phi_a = obj.phi_a(logintMu);
% % Shear rate corresponding to steady shear structure 
% A =  - obj.par.b_0 * obj.cnst.phi_p ...
%     * ((phi_a/obj.cnst.phi_p)^((2)/(obj.par.d_f - 3)) ...
%     -  obj.par.m_p^(1/obj.par.d_f)*(phi_a/obj.cnst.phi_p)^((3)/(obj.par.d_f - 3)));
% 
% B = 4 * cutOff(obj, logintMu) * obj.par.alfa ...
%     * (obj.cnst.phi_p^2/pi) ...
%     * (phi_a/obj.cnst.phi_p)^((obj.par.d_f)/(obj.par.d_f - 3));
% 
% phi_a = obj.phi_a(logintMu);
% eta = obj.etaTrimodal(logintMu);
% 
% C = 2 * cutOff(obj, logintMu) * obj.cnst.k_b * obj.cnst.T * obj.cnst.phi_p^2 ...
%    ./( 2 * obj.cnst.mu_s * eta * obj.par.W * pi * obj.cnst.a_p^3) ...
%    * (phi_a/obj.cnst.phi_p).^((obj.par.d_f + 3)/(obj.par.d_f - 3));
% 
% structure_SR = abs((B + sqrt(B.^2-4.*A.*C))./(2.*A));
% x = structure_SR;


%%

%Model Parameters
a_p   = obj.cnst.a_p;    % primary particle radius
phi_p = obj.cnst.phi_p;  % particle volume fraction
k_b   = obj.cnst.k_b;    % Boltzmann constant
T     = obj.cnst.T;      % Temperature
mu_s  = obj.cnst.mu_s;   % Medium viscosity (Pa-s)

df    = obj.par.d_f;     % Fractal dimension     
alfa  = obj.par.alfa;    % Collision efficiency
W     = obj.par.W;       % Fuch's stability ratio
b_o   = obj.par.b_0;     % Breakage constant
m_p   = obj.par.m_p;     % Number of primary particles in a primary cluster


intMu=exp(logintMu); 


fra_moment = @(j, c) obj.fra_moment(j, c);
theta = @(x) obj.theta(x);
shearcoe=phi_p*alfa*obj.cutOff(logintMu)/(2*pi);

browncoe=phi_p*k_b*T*obj.cutOff(logintMu)...
    /(4*pi*(a_p^3)*mu_s*obj.etaTrimodal(logintMu)*W);

breakcoe= b_o;

c=obj.MOMIC(logintMu);

 % 0th moment 
    SA = -2*shearcoe*...
        (intMu(1)*fra_moment(3/df,c)+...
        3*fra_moment(2/df,c)*fra_moment(1/df,c));
    
    BR = -2*browncoe*...
        (intMu(1)^2+fra_moment(1/df,c)*fra_moment(-1/df,c));
    
    SB = breakcoe*(theta(0)-1)*...
        (fra_moment(1/df,c)-m_p^(1/df)*intMu(1));
    
structure_SR = abs((-SA + sqrt(SA.^2 - 4*SB*BR))./2/SB);
x = structure_SR;
end