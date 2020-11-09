function out = UDLAOS(obj, gamma_0, omega, time)
% UD-LAOS stress calculation

% Sinusoidal shear rate
%shear_rate = @(t) gamma_0*omega + gamma_0*omega*cos(omega*t);
shear_rate =@(t) gamma_0*omega;

% Scaled number density
loginitialMu=obj.InitialCondition();
% Setting tolerance for ODEs
odeopts = odeset('RelTol',1e-3,'AbsTol', 1e-6, 'Stats','on');

% Solving for transient state at specified times
fun = @(t, X) [obj.momicDerivative9(t, X(1:9)', shear_rate(t)); 
               obj.elasticStrain(t, X(end), shear_rate(t), X(1:9)')];
           
out.sol = ode15s(fun, time, [loginitialMu, obj.gamma_lin], odeopts);
% 
% try
%     % Obtaining solution variables
%     [sol_Eval, derivatives] = deval(sol, tDim);
%     out.phi_a = obj.cnst.phi_p.*(sol_Eval(1,:)'./obj.par.m_p).^(1-3/obj.par.d_f);
%     out.gamma_e =  min( obj.gamma_lin, sol_Eval(2,:))';
%     out.gamma_e_dot = derivatives(2,:)';
% 
%     [~, msgid] = lastwarn;
%     if strcmp(msgid, 'MATLAB:illConditionedMatrix')
%         error('Matrix is singular, close to singular or badly scaled. Results may be inaccurate')
%     end
%     lastwarn('')
% catch
%     warning('Problem evaluating ODE.  Assigning a value of 0.');
%     out.phi_a = 10^6*ones(size(tDim));
%     out.gamma_e = 10^6*ones(size(tDim));
%     lastwarn('')
%     out.EXITFLAG = -10;
% end
% 
%     % Computing overall macroscopic stress
%     out.stress = total_stress(obj, phi_a, gamma_e, shear_rate(time));
%     end