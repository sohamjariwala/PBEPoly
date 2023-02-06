function [L, W] = quadrature_solve(obj,logintMu)
% QUADRATURE_SOLVE Convert moments to abscissa and weights.
% The expression used in this function is based on the product-difference
% algorithm continued fraction approximation of moments due to Gordon (1968).

mu = exp([logintMu(1) 0 logintMu(2:5)])./exp(logintMu(1));

% Elements of the Jacobi matrix
a1 = mu(2);
b1 = sqrt(mu(3) - mu(2)^2);
a2 = (mu(2)^3 - 2*mu(2)*mu(3) + mu(4))/(mu(3) - mu(2)^2);
b2 = sqrt((-mu(3)^3 + 2*mu(2)*mu(3)*mu(4) - mu(4)^2 ...
      - mu(2)^2*mu(5) + mu(3)*mu(5))/(mu(3) - mu(2)^2)^2);

a3 = (mu(4)^3 - 2*mu(3)*mu(4)*mu(5)...
     - mu(2)*mu(3)*(mu(3)^3 + 2*mu(4)^2 - 2*mu(3)*mu(5)) - mu(2)^3*(mu(4)^2 + 2*mu(3)*mu(5))...
     + mu(2)^4*mu(6) + mu(3)^2*mu(6) + mu(2)^2*(3*mu(3)^2*mu(4) + 2*mu(4)*mu(5) -2*mu(3)*mu(6))) ...
     /((mu(2)^2 - mu(3))*(mu(3)^3 + mu(4)^2 + mu(2)^2*mu(5) - mu(3)*(2*mu(2)*mu(4) + mu(5))));

% Tridiagonal Jacobi matrix
J = [a1, b1,  0;
     b1, a2, b2;
      0, b2, a3];

% Eigenvalues of the Jacobi matrix
[V,D] = eig(J);

% Weights and Absicissa
L = [D(1,1); D(2,2); D(3,3)].^(1/obj.par.d_f)*2*obj.cnst.a_p;
W = mu(1)*[V(1,1); V(1,2); V(1,3)].^2;

end