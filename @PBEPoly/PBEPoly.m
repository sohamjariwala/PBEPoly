classdef PBEPoly
% A rheological constitutive model thixo-elasto-viscoplastic materials
% via population balance modeling.

%% Parameters
    properties 
       % Parameters to be determined from fitting  
       par = struct(...
            'W', [], ...                % Fuch's stability ratio
            'alfa', [], ...             % Collision efficiency
            'b_0' , [], ...             % Breakage constant
            'd_f' , [], ...             % Fractal dimension
            'porosity', [], ...         % R_h/R_a
            'm_p', [],...               % #(particles) in a primary cluster
            'kh',[],...                 % Backstress modulus parameter
            'p',[]...                   % Exponent for relaxation time
            );
  
       % Constants that will remain fix throughout the computation.
       cnst = struct(... 
            'phi_p'    , [], ...        % Particle volume fraction
            'a_p'      , [], ...        % Primary particle radius
            'mu_s'     , [], ...        % Suspension viscosity (Pa-s)
            'sigma_y0' , [], ...        % Yield stress (Pa)
            'phi_max'  , 0.64, ...      % Maximum packing fraction
            'k_b' , 1.38064852e-23, ... % Boltzmann constant
            'T' , 298.15 ...            % Temperature
        );
       % Inverse of the matrix that scales the logarithm of moments
        lnvA;

       % Minimum allowed agglomerate volume fraction
        phi_pc;
    end
    
%% Functions describing rheological variables
methods
    function obj = PBEPoly()
        % Initialize fitting parameters
        obj.par.W = 0.008998930690468;
        obj.par.alfa = 0.665487794544652;
        obj.par.b_0 = 0.003742494709682;
        obj.par.d_f = 2.11;
        obj.par.porosity = 0.92;
        obj.par.m_p = 468;
        obj.par.kh = 0.7;
        obj.par.p = 3;

        % Initialize constants
        obj.cnst.phi_p = 0.03;
        obj.cnst.a_p = 8e-9;
        obj.cnst.mu_s = 0.41;
        obj.cnst.sigma_y0 = 11;
        obj.cnst.G_0 = 450;
        obj.cnst.phi_max = 0.64;
        obj.cnst.q = 1;
        obj.cnst.k_b = 1.3806e-23;
        obj.cnst.T = 298.15;
        obj = inverseOfA(obj);
    end

    function x = get.phi_pc(obj)
        % GET.PHI_PC Minimum allowable volume fraction of aggregates
        x = obj.cnst.phi_p*obj.par.m_p^(3/obj.par.d_f-1);
    end

    function x = phi_h(obj, phi)
        % PHI_H Hydrodynamic volume fraction
        x = obj.par.porosity^3.*phi;
    end

    function x = phi_a(obj, logintMu)
        % PHI_A Volume fraction computed from the moments of the distribution
        c = obj.MOMIC(logintMu);
        x = obj.cnst.phi_p*obj.fra_moment(3/obj.par.d_f,c);
    end

    function x = cutOff(obj, logintMu)
        % CUTOFF Cutoff function to model dynamic arrest
        phi_a = obj.phi_a(logintMu);

        x = abs((obj.phi_max(logintMu) - phi_a)...
            ./ (obj.phi_max(logintMu) - obj.phi_pc))^(2/(3-obj.par.d_f));
    end

    function x = phi_max(~, logintMu)
        % PHI_MAX Maximum packing fraction as a function of moments (polydisperse)
        Mu = exp(logintMu); Mu = [Mu(1) 1 Mu(2:end)];
        var_sigma = sqrt(abs(reallog( ...
            Mu(4)*Mu(2)/Mu(3)^2 ...
            )));

        x =  1-0.57*exp(-var_sigma) + 0.2135*exp(-0.57*var_sigma/0.2135) ...
            + 0.0019*(cos(2*pi*(1-exp(-0.75*var_sigma.^0.7 - 0.0025*var_sigma.^4)))-1);
    end

    function x = eta(obj, phi, logintMu)
        % ETA Suspension viscosity based on Krieger-Dougherty equation.
        phi_max = obj.phi_max(logintMu);
        x = abs( 1 - obj.phi_h(phi) / phi_max ).^(-2.5*phi_max);
    end

    function x = sigma_y(obj, logintMu)
        % SIGMA_Y Elastic modulus based on Shih et al. strong link regime
        d_b = 1; % backbone fractal dimension
        x = obj.cnst.sigma_y0*abs((obj.phi_a(logintMu) - obj.phi_pc)...
            ./(obj.phi_max(logintMu) - obj.phi_pc)).^(2/(3-obj.par.d_f));
    end

    function x = gamma_dot_p(obj,sigma,A,logintMu,~)
        % GAMMA_DOT_P Plastic deformation rate

        x = abs(obj.sigma_eff(sigma,logintMu,A)) ...
            /(obj.sigma_y(logintMu)/obj.structure_shear_rate(logintMu) ...
            +obj.cnst.mu_s*obj.etaTrimodal(logintMu));

    end

    function x = viscous_stress(obj, logintMu, shear_rate)
        % Viscous stress
        x = obj.cnst.mu_s.*obj.etaTrimodal(logintMu).*shear_rate;
    end

    function x = total_stress_SS(obj, logintMu, shear_rate)
        % Total stress
        x = (obj.par.kh/obj.cnst.q + obj.sigma_y(logintMu))*sign(shear_rate) ...
            + obj.viscous_stress(logintMu, shear_rate);
    end

    function x = tau(obj, logintMu)
        % TAU Viscoelastic relaxation time
        x = abs(obj.cnst.sigma_y0/obj.sigma_y(logintMu))^-obj.par.p ...
        *abs(obj.sigma_y(logintMu)/obj.cnst.G_0...
            ./obj.structure_shear_rate(logintMu));
    end

    function x = Adot(obj, ~, A, logintMu, sigma,shear_rate)
        % ADOT Derivative of kinematic hardening parameter
        gamma_dot_p = sign(obj.sigma_eff(sigma, logintMu, A)) ...
        *obj.gamma_dot_p(sigma,A,logintMu,shear_rate);
        x = gamma_dot_p - obj.cnst.q*A*abs(gamma_dot_p);
    end

    function x = sigma_eff(obj, sigma, logintMu, A)
        % SIGMA_EFF Effective stress after kinematic hardening
        x = sigma - obj.par.kh*(obj.sigma_y(logintMu)/obj.cnst.sigma_y0)*A;
    end

    %% Constitutive models

    % Viscosity for a polydisperse suspension.
    x = etaTrimodal(obj, logintMu)

    % Shear rate corresponding to steady shear structure
    x = structure_shear_rate(obj, logintMu)

    % Shear stress DE corresponding to the Saramito viscoelastic model
    dX = shearStressDE(obj, t, sigma, gamma_e, shearRate, logintMu)

    %% Stress responses
    out = steadyShear(obj, shearRate, initialConditions)

    out = steadyShearODE(obj, shearRate, initialConditions)

    out = stepShear(obj, initialShearRate, finalShearRate, time, init)

    out = UDLAOS(obj, gamma_0, omega, time, init)

    out = stressResponse(obj, initialShearRate, shearRate, time, initialConditions)
    
    out = shearReversal(obj, shear_rate, time, initialConditions)

    out = shearRateJump(obj, init_rate, delta_t, min_rate, final_rate, time, initialConditions)

    out = LAOS(obj, gamma_0, omega, time, initialConditions)

    %% Microstructure
    % MOMIC
    function x = MOMIC(obj, logintMu)
        %Using MOMIC Method to find the coefficients from integer moments

        x=obj.lnvA*(logintMu');

    end

    function obj = inverseOfA(obj)
        % Create the A matrix for
        K=5;
        A=zeros(K,K);
        j=[0,2:K];
        k=j;
        for a=1:K
            for b=1:K
                A(a,b)=j(a)^k(b)-j(a);
            end
        end
        obj.lnvA=inv(A);
    end

    function x = fra_moment(~,j,c)
        % Fractional moments
        K=length(c);
        lnintMu_k=c(1)*(1-j);
        for i=2:K
            lnintMu_k=lnintMu_k+c(i)*(j^i-j);
        end
        x=exp(lnintMu_k);
    end

    function x=theta(~, k)
        % Moment of the daughter distribution
        r=1;
        p=2; % #(Daughter fragments)
        q=r/(p-1);
        x=p*gamma(q+k)/gamma(q)*gamma(q+r)/gamma(q+r+k);
    end

    function x=InitialCondition(obj)
        % Default initial guess for the 5 moments
        mn=obj.par.m_p;
        %assuming initially monodisperse
        gam=[1,1,1,1,1,1];
        K=length(gam);
        initialMu=zeros(size(gam));
        for k=1:K
            initialMu(k)=gam(k)*(mn^(k-2));
        end
        initialMu(2)=[];
        x=reallog(initialMu);

        x = [-6.221327387190237 6.411211053369301 13.126103089903342 ...
            20.195599812988782 27.608183731624543];
    end

    % Microstructure evolution through population balances
    f = momicDerivative5(obj, t, shearRate,logintMu);

    %% Plot functions
    function flowcurve(obj,shear_rate)
        % FLOWCURVE Plot steady state flow curve
        if nargin < 2
            shear_rate = logspace(-2,2,20);
        end

        for i = length(shear_rate):-1:1
            if i == length(shear_rate)
                out = obj.steadyShear(shear_rate(i));
                stress(i) = out.stress;
                logintMu(i,:) = out.logintMu;
                elastic_stress(i) = out.sigma_y;
            else
                out = obj.steadyShear(shear_rate(i), out);
                stress(i) = out.stress;
                logintMu(i,:) = out.logintMu;
                elastic_stress(i) = out.sigma_y;
            end
        end

        figure('Name','Steady state flow curve','NumberTitle','off')
        loglog(shear_rate, stress, ...
            shear_rate, elastic_stress, ...
            shear_rate, stress-elastic_stress, ...
            'LineWidth', 2);
        xlabel('$\dot\gamma \ [\mathrm{s}^{-1}]$','Interpreter','latex');
        ylabel('$\sigma_{yx} \ [\mathrm{Pa}]$','Interpreter','latex');
        legend('Shear stress','Elastic stress','Viscous stress', ...
            'Location','best')
        set(gca,'FontSize',20,'LineWidth',2)
        axis([-inf inf 0.1 inf])
    end

    function distribution(obj, caxis, logintMu)
        %% Distribution generation function for MOMIC
        % This function generates lognormal distribution for provided physical
        % moments in log scale.
        % INPUT: Log scale moments
        % OUTPUT: distribution function and plot

        c=obj.MOMIC(logintMu);

        gamma2 = exp(logintMu(:,1)).*exp(logintMu(:,2));

        a_p = obj.cnst.a_p;
        d_f = obj.par.d_f;
        N = 1000; % Number of points

        figure('Name','Probability density of agglomerate sizes')
        c_map = colormap(copper(length(caxis)));
        lnsigmag = zeros(size(caxis));
        mulnX = zeros(size(caxis));
        for i = 1:length(caxis)
            lnsigmag(i) = sqrt(log(gamma2(i)));
            mulnX(i) = 1/sqrt(gamma2(i));
            average(i) = exp(mulnX(i) + lnsigmag(i)^2/2);
            a = linspace(-8,8,N);
            b =  log(10)*(exp(a).^(1/d_f))...
                .*(1/(sqrt(2*pi)*lnsigmag(i))...
                *exp(-0.5*((a - log(mulnX(i)))/(lnsigmag(i))).^2));
            a = 10^6*(exp(a)*obj.fra_moment(1+1/d_f,c(:,i))*2*a_p);
            semilogx(a', b', 'Color', c_map(i,:), 'LineWidth',2);
            hold on
        end

        xlabel('$\mathrm{Aggregate\ size\ (diameter)} (\mu m)$','Interpreter','latex');
        ylabel('$\ln(10) x^{1+1/d_f} f(x)$','Interpreter','latex');
        axis([-inf, inf, 0, inf])
        axis square
        clim([caxis(1) caxis(end)])
        set(gca,'FontSize',14,'FontWeight','bold','linewidth',2);
        colorbar, grid on
    end

    %% Event and helper functions
    function [values,isterminal,direction] = odeEvent(obj, t,X,tstart)
        %  ODEEVENT Event function to stop ODE calculation after 3 s
        values(1) = t;
        %  Don't let integration go for more than 3 seconds.
        values(2) = toc(tstart) < 10;
        isterminal = true(size(values));
        direction = zeros(size(values));
    end
end
end