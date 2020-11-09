    function out = stressResponse(obj, initialShearRate, shearRate, time)
    % Stress response calculation

    % shear rate as a function of time
    shear_rate = @(t) interp1(time, shearRate, t, 'linear');

    [stress, phi_a_init, gamma_e_init, EXITFLAG] = steadyShear(obj, initialShearRate);

    % Scaled number density
    nu_init = (phi_a_init/obj.cnst.phi_p)^(obj.par.d_f/(obj.par.d_f-3));

    % Setting tolerance for ODEs
    %odeopts = odeset('RelTol',1e-3,'AbsTol', [1e-6 1e-6], 'Stats','on');

    % Dimensionless time span
%     tDim = time./obj.par.b_0;
    tDim = time;

    % Solving for transient state at specified times
    fun = @(t, X) mwasameModel(obj, t, X, shear_rate(t));
    sol = ode23tb(fun, tDim, [nu_init, obj.gamma_lin]);

    try
        % Obtaining solution variables
        [sol_Eval derivatives] = deval(sol, tDim);
        out.phi_a = obj.cnst.phi_p.*(sol_Eval(1,:)'./obj.par.m_p).^(1-3/obj.par.d_f);
        out.gamma_e =  min( obj.gamma_lin, sol_Eval(2,:))';
        out.gamma_e_dot = derivatives(2,:)';
        [msgstr, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:illConditionedMatrix')
            error('Matrix is singular, close to singular or badly scaled. Results may be inaccurate')
        end
        warning('')

    catch
        warning('Problem evaluating ODE.  Assigning a value of 0.');
        out.phi_a = 10^6*ones(size(tDim));
        out.gamma_e = 10^6*ones(size(tDim));
        warning('')
    end
     % Computing overall macroscopic stress
    out.stress = total_stress(obj, phi_a, gamma_e, shear_rate(time));    
    end
    