classdef PBEPoly
% A rheological constitutive model for carbon black via population
% balances.
    
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
            'kh',[]...                  % Backstress modulus parameter
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
        lnvA;
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
            obj.par.kh = 0;

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
            obj.phi_pc =  obj.calculate_phi_pc();
        end

        function x = calculate_phi_pc(obj)
        % Minimum allowable volume fraction of aggregates
            x = obj.cnst.phi_p*obj.par.m_p^(3/obj.par.d_f-1);
        end

        function x = phi_h(obj, phi)
        % Hydrodynamic volume fraction
            x = obj.par.porosity^3.*phi;
        end
        
        function x = phi_a(obj, logintMu)
        % Volume fraction computed from the moments of the distribution
            c = obj.MOMIC(logintMu);
            x = obj.cnst.phi_p*obj.fra_moment(3/obj.par.d_f,c);
        end
                
        function x = cutOff(obj, logintMu)
        % Cutoff function to model dynamic arrest
            phi_a = obj.phi_a(logintMu);

            x = abs((obj.phi_max(logintMu) - phi_a)...
                ./ (obj.phi_max(logintMu) - obj.phi_pc))^(2/(3-obj.par.d_f));
        end
        
        function x = phi_max(~, logintMu)
        % Maximum packing fraction as a function of moments (polydisperse)
            Mu = exp(logintMu); Mu = [Mu(1) 1 Mu(2:end)];
            var_sigma = sqrt(abs(reallog( ...
                Mu(4)*Mu(2)/Mu(3)^2 ...
                )));
   
            x =  1-0.57*exp(-var_sigma) + 0.2135*exp(-0.57*var_sigma/0.2135) ...
            + 0.0019*(cos(2*pi*(1-exp(-0.75*var_sigma.^0.7 - 0.0025*var_sigma.^4)))-1);
        end
        
        function x = eta(obj, phi, logintMu)
        % Suspension viscosity based on Krieger-Dougherty equation.
            phi_max = obj.phi_max(logintMu);
            x = abs( 1 - obj.phi_h(phi) / phi_max ).^(-2.5*phi_max);
        end
        
        % Viscosity for a polydisperse suspension.
        x = etaTrimodal(obj, logintMu)
        
        function x = sigma_y(obj, logintMu)
        % Elastic modulus based on Shih et al. strong link regime
        d_b = 1; % backbone fractal dimension
        x = obj.cnst.sigma_y0*abs((obj.phi_a(logintMu) - obj.phi_pc)...
         ./(obj.phi_max(logintMu) - obj.phi_pc)).^(2/(3-obj.par.d_f));
        end
        
        function x = gamma_dot_p(obj, sigma, A, logintMu)
        % Plastic deformation rate
            if abs(obj.sigma_eff(sigma, A)) < obj.sigma_y(logintMu)
                x = 0;

            elseif  abs(obj.sigma_eff(sigma, A)) >= obj.sigma_y(logintMu)
                x = (abs(obj.sigma_eff(sigma, A)) - obj.sigma_y(logintMu))...
                    * sign(obj.sigma_eff(sigma, A))/obj.etaTrimodal(logintMu);
            end
        end
        
        function x = viscous_stress(obj, logintMu, shearRate)
        % Viscous stress
            x = obj.cnst.mu_s.*obj.etaTrimodal(logintMu).*shearRate;
        end
        
        function x = total_stress_SS(obj, logintMu, shearRate)
        % Total stress
           x =  obj.sigma_y(logintMu) ...
              + obj.viscous_stress(logintMu, shearRate);
        end
        
        function x = tau(obj, logintMu)
        % Relaxation time      
            x = abs(obj.cnst.sigma_y0/obj.cnst.G_0...
                ./obj.structure_shear_rate(logintMu));
        end
     
        %% Stress responses
        out = steadyShear(obj, shearRate, initialConditions)
        
        out = steadyShearODE(obj, shearRate, initialConditions)
        
        out = stepShear(obj, initialShearRate, finalShearRate, time, init)
        
        out = UDLAOS(obj, gamma_0, omega, time, init)
        
        out = stressResponse(obj, initialShearRate, shearRate, time)

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
         
         %% Constitutive models
         % Shear rate corresponding to steady shear structure 
         x = structure_shear_rate(obj, shearRate, logintMu)

         % Shear stress DE corresponding to the Saramito viscoelastic model
         dX = shearStressDE(obj, t, sigma, gamma_e, shearRate, logintMu)

         function x = Adot(obj, ~, A, logintMu, sigma)
             % ADOT Derivative of kinematic hardening parameter
             x = obj.gamma_dot_p(sigma,A,logintMu) ...
                 - obj.cnst.q*A*abs(obj.gamma_dot_p(sigma,A,logintMu));
         end

         function x = sigma_eff(obj, sigma, A)
             % SIGMA_EFF Effective stress after kinematic hardening
             x = sigma - obj.par.kh*A;
         end
    end
end