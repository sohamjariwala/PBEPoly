function eta = etaTrimodal(obj, logintMu)
% ETATRIMODAL Calculate viscosity for polydisperse colloidal systems
% based on the methodology described in Mwasame, Wagner, Beris, "Modeling the 
% effects of polydispersity on the viscosity of noncolloidal hard sphere
% suspensions"
 
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

    beta_star = beta(d,DD, phi_d, phi_DD);
    alfa_star = (beta_star*phi_D*D + (1 - beta_star)*phi_d*d)/(beta_star*phi_D + (1 - beta_star)*phi_d);

    beta1 = beta(d, alfa_star, phi_d, phi_D + phi_DD);
    beta2 = beta(alfa_star, DD, phi_d + phi_D, phi_DD);

    fTri = fu(beta2*(beta1*phi_d+phi_D)+phi_DD) ...
         + fu(beta1*phi_d+phi_D)*(1-beta2) ...
         + fu(phi_d)*(1-beta1);
    
     eta = real(exp(fTri));

end