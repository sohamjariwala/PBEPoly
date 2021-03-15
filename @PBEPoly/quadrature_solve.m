function [L, W] = quadrature_solve(obj, logintmu)
% Function to calculate the the equivalent ternary system using the
% following decomposition:
%        3
% m_k = sum omega (V_i)^k
%       i=1
% where V_i is the volume scaled by the number mean.
% -------------------------------------------------------------------------
% Log-normal test (Mwasame et al. (2016a))
% mu = 4.58; sigma = 0.36;
% for i = 0:5
%     m(i+1) = exp(i*mu + i^2*sigma^2);
% end

% Inserting log(mu_1)
logintmu = [logintmu(1) 0 logintmu(2:end)];

mu = exp(logintmu);

% Normalizing the scaled moments
gamma = zeros(size(mu));
for i = 1:length(mu)
   gamma(i) = mu(i)/mu(1)^(2-i)/mu(2)^(i-1);
end

solveThis = @(x)  jQuadsolve(x, gamma);

% Nonlinear solver
options = optimset('TolFun', eps, 'TolX', eps, 'Display', 'off');
optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true);

[x_sol,~,~,~] = ...                      %[x_sol,RESNORM,RESIDUAL,EXITFLAG]
lsqnonlin(solveThis, ...
  [32.3027; ...
    0.2894; ...
  530.8157; ...
    0.0238; ...
    0.9737; ...
    0.0000], ...  % Init value
       zeros(6,1), ...                            % Lower bound
       [inf inf inf 1 1 1]', ...                  % Upper bound
       options);

  R = mu(1)^(-1/obj.par.d_f)*obj.cnst.a_p;
  radius = x_sol(1:3).^(1/3)*R;
  weights = x_sol(4:end);

  A = sortrows([radius weights], 1);
  L = A(:,1);
  W = A(:,2);
  
function [F2, Jacobian] = jQuadsolve(x, gamma)
F2 =     ([    1,      1,      1;
              x(1)^1, x(2)^1, x(3)^1;
              x(1)^2, x(2)^2, x(3)^2;
              x(1)^3, x(2)^3, x(3)^3;
              x(1)^4, x(2)^4, x(3)^4;
              x(1)^5, x(2)^5, x(3)^5] * [x(4); x(5); x(6)])./gamma(1:6)' - 1;

    if nargout > 1
        Jacobian = zeros(6,6);    
        for i = 1:6
                Jacobian(i,1:3) = (i-1)*[x(4) x(5) x(6)]/gamma(i);
                Jacobian (i, 4:6) = [x(1) x(2) x(3)].^(i-1)/gamma(i);
        end
    end
end

  
end