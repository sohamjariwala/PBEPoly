function eta = etaTrimodal(obj, logintMu)
% This function calculates the viscosity for polydisperse systems based on
% the methodology described in Mwasame, Wagner, Beris, "Modeling the 
% effects of polydispersity on the viscosity of noncolloidal hard sphere
% suspensions"

try 
    [L, W] = quadrature_solve(obj, logintMu);

    d = L(1); D = L(2); DD = L(3);

    phi = obj.phi_a(logintMu);

    phi_d  = W(1)/sum(W) * phi;
    phi_D  = W(2)/sum(W) * phi;
    phi_DD = W(3)/sum(W) * phi;

    % Equation (36) in Mwasame et al.
    beta = @(d, D, phi_d, phi_D) ...
        ((d/D)^0.18)^(1-(0.54*phi_d/(0.54*phi_d + phi_D)));

    % Equation (2) the obj.eta function calls the Krieger-Dougherty
    fu = @(x) log(obj.eta(x, logintMu));

    beta1 = beta(d, D, phi_d, phi_D);
    beta2 = beta(D, DD, phi_D, phi_DD);

    fTri = fu(beta2*(beta1*phi_d+phi_D)+phi_DD) ...
         + fu(beta1*phi_d+phi_D)*(1-beta2) ...
         + fu(phi_d)*(1-beta1);
    
     eta = exp(fTri);
catch
    eta = 10^6; 
end
end