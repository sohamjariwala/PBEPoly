function out = LAOS(obj, gamma_0, omega, time, initialConditions)
% UD-LAOS stress calculation
% Solve for stress and structural variables for step change in shear rate

if nargin > 4
    init.logintMu = initialConditions.logintMu;
    init.stress = initialConditions.stress;
    init.A = initialConditions.A;
else
    init = obj.steadyShearODE(initialShearRate);
end

% Number of moments to solve for
numMoments = 5;

% Step-up and Step-down
shear_rate = @(t) gamma_0*omega*cos(omega*t);

%% Solving the system of ODEs
% Setting tolerance for ODEs
tstart = tic;
odeopts = odeset('RelTol',1e-5,'AbsTol',1e-6,'Stats','off','Events', ...
    @(t,X) obj.odeEvent(t,X,tstart));

% Simultaneous ODEs to be solved till stationary state is attained
fun = @(t, X) [obj.momicDerivative5(t, X(1:numMoments)', obj.gamma_dot_p(X(end), X(end-1), X(1:numMoments)',shear_rate(t)));
    obj.Adot(t,X(end-1), X(1:numMoments)',X(end),shear_rate(t));
    obj.shearStressDE(t, X(end), shear_rate(t), X(1:numMoments)',X(end-1))];

% try
    %% Running the solution
    out.sol = ode15s(fun, ...
        [0, time(end)], ...
        [init.logintMu, init.A, init.stress],...
        odeopts);

    % Obtaining solution variables
    if length(time) <= 2
        time = out.sol.x;
    end

    [out.sol_Eval, out.derivatives] = deval(out.sol, time);
    out.time = time;
    out.logintMu = out.sol_Eval(1:numMoments,:)';
    out.A = out.sol_Eval(end-1,:);
    out.stress = out.sol_Eval(end,:);
    for i = 1:length(time)
        out.phi_a(i) = obj.phi_a(out.logintMu(i,:));
        out.sigma_y(i) = obj.sigma_y(out.logintMu(i,:));
    end
    out.EXITFLAG = 1;
    out.strain = gamma_0*sin(omega*time);
    out.shear_rate = gamma_0*omega*cos(omega*time);

% catch
%     out.time = time;
%     out.logintMu = 10^6*ones(length(time),1)*init.logintMu;
%     out.A = 10^6*ones(1,length(time));
%     out.stress = 10^6*ones(1,length(time));
%     out.phi_a = 10^6*ones(1,length(time));
%     out.sigma_y = 10^6*ones(1,length(time));
%     out.EXITFLAG = -10;
%     return
% end
end