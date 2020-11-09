function out = stepShear(obj, initialShearRate, finalShearRate, time)
    % Step-up and Step-down
    shear_rate = @(t) finalShearRate;

    % Value at initial steady state
    [~, phi_a_init, gamma_e_init, out.EXITFLAG] = steadyShear(obj, initialShearRate);

    mu_init = (phi_a_init/obj.cnst.phi_p)^(obj.par.d_f/(obj.par.d_f-3));

    % Setting tolerance for ODEs
    odeopts = odeset('RelTol',1e-8,'AbsTol', [1e-9 1e-9], 'Stats','on');

    % Solving for transient state at specified times
    fun = @(t, X) mwasameModel(obj, t, X, shear_rate(t));
    sol = ode23tb(fun, time, [mu_init, gamma_e_init], odeopts);

    % Obtaining solution variables
    sol_Eval = deval(sol, time);
    out.phi_a = obj.cnst.phi_p.*(sol_Eval(1,:)').^(1-3/obj.par.d_f);
    out.gamma_e = min( obj.gamma_lin, sol_Eval(2,:))';

    % Computing overall macroscopic stress
    out.stress = total_stress(obj, phi_a, gamma_e, shear_rate(time));
end
