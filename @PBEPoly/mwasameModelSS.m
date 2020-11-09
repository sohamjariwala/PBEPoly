function [d_phi_a] = mwasameModelSS(obj, phi_a, shear_rate)
% Microstructure at steady state through population balances   
    % Brownian aggregation term
    brownianAggregation = ...
        2 * cutOff(obj, phi_a) * obj.cnst.k_b * obj.cnst.T * obj.cnst.phi_p^2 ...
       ./( 2 * obj.cnst.mu_s * eta(obj, phi_a) * obj.par.W * pi * obj.cnst.a_p^3) ...
       .* (3-obj.par.d_f)/obj.par.d_f ...
       .* (phi_a/obj.cnst.phi_p)^((2*obj.par.d_f + 3)/(obj.par.d_f - 3));

    % Shear aggregation term
    shearAggregation = ...
        4 * cutOff(obj, phi_a) * obj.par.alfa ...
        .* (3-obj.par.d_f)/obj.par.d_f ...
        .* (obj.cnst.phi_p^2/pi) ...
        .* (phi_a/obj.cnst.phi_p).^((2*obj.par.d_f)/(obj.par.d_f - 3)) ...
        .* abs(shear_rate);


    % Shear breakage term
    shearBreakage = ...
        - obj.par.b_0 * obj.cnst.phi_p ...
        * (3-obj.par.d_f)/obj.par.d_f ...
        * ((phi_a/obj.cnst.phi_p)^((obj.par.d_f + 2)/(obj.par.d_f - 3)) ...
        - (obj.par.m_p)^(1/obj.par.d_f)*(phi_a/obj.cnst.phi_p)^((obj.par.d_f + 3)/(obj.par.d_f - 3))) ...
        * abs(shear_rate).^2;

        d_phi_a = 10;
        if (shearBreakage~=0)
         d_phi_a = (brownianAggregation + ...
              shearAggregation + ...
              shearBreakage)/(shearBreakage);
        end
end