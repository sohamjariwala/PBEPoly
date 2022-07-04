classdef PBEPoly
% A rheological constitutive model for carbon black via population
% balances.
    
%% Parameters    
    properties 
       % Parameters to be determined from fitting  
       par = struct(...
            'W', [], ...                % Fuch's stability ration
            'alfa', [], ...             % Collision efficiency
            'b_0' , [], ...             % Breakage constant
            'd_f' , [], ...             % Fractal dimension
            'porosity', [], ...         % R_h/R_a
            'm_p', []...                % #(particles) in a primary cluster
            );
  
       % Constants that will remain fix throughout the computation.
       cnst = struct(... 
            'phi_p'    , [], ...        % particle volume fraction
            'a_p'      , [], ...        % primary particle radius
            'mu_s'     , [], ...        % Suspension viscosity (Pa-s)
            'G_0'      , [], ...        % Equilibrium modulus (Pa)
            'sigma_y0' , [], ...        % Yield stress (Pa)
            'phi_max'  , 0.64, ...      % Maximum packing fraction
            'k_b' , 1.38064852e-23, ... % Boltzmann constant
            'T' , 298.15 ...            % Temperature
        );
        lnvA;
        phi_pc;
        gamma_lin;
    end
%% Functions describing rheological variables
    methods
        function obj = PBEPoly(obj)
          obj.par.W = 0.008998930690468;
          obj.par.alfa = 0.665487794544652;
          obj.par.b_0 = 0.003742494709682;
          obj.par.d_f = 2.11;
          obj.par.porosity = 0.92;
          obj.par.m_p = 468;
          
          obj.cnst.phi_p = 0.03;
          obj.cnst.a_p = 8e-9;
          obj.cnst.mu_s = 0.41;
          obj.cnst.G_0 = 450;
          obj.cnst.sigma_y0 = 11;
          obj.cnst.phi_max = 0.64;
          obj.cnst.k_b = 1.3806e-23;
          obj.cnst.T = 298.15;
          obj = inverseOfA(obj);
        end
                
        function x = get.gamma_lin(obj)
        % Limit of linearity of elastic strain
            x = obj.cnst.sigma_y0/obj.cnst.G_0;
        end
        
        function x = get.phi_pc(obj)
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
        % Hyperbolic cutoff function to model dynamic arrest
            phi_a = obj.phi_a(logintMu);

       x = ((obj.phi_max(logintMu) - phi_a)...
         ./ (obj.phi_max(logintMu) - obj.phi_pc))^(2/(3-obj.par.d_f));
        
        end
        
        function x = phi_max(~, logintMu)
        % Maximum packing fraction as a function of moments (polydisperse)
            Mu = exp(logintMu); Mu = [Mu(1) 1 Mu(2:end)];
            sigma = sqrt(log( ...
                Mu(4)*Mu(2)/Mu(3)^2 ...
                ));
   
            x =  1-0.57*exp(-sigma) + 0.2135*exp(-0.57*sigma/0.2135) ...
            + 0.0019*(cos(2*pi*(1-exp(-0.75*sigma.^0.7 - 0.0025*sigma.^4)))-1);
        end
        
        function x = eta(obj, phi, logintMu)
        % Suspension viscosity based on Krieger-Dougherty equation.
            phi_max = obj.phi_max(logintMu);
            x = ( 1 - obj.phi_h(phi) / phi_max ).^(-2.5*phi_max);
        end
        
        % Viscosity for a polydisperse suspension.
        x = etaTrimodal(obj, logintMu)
        
        function x = G(obj, logintMu)
        % Elastic modulus based on Shih et al. strong link regime
        d_b = 1; % backbone fractal dimension
        x = obj.cnst.G_0.*((obj.phi_a(logintMu) - obj.phi_pc)...
         ./(obj.phi_max(logintMu) - obj.phi_pc)).^((2)/(3-obj.par.d_f));
        end
        
        function x = elastic_stress(obj, logintMu, gamma_e)
        % Elastic stress
            x = obj.G(logintMu).*gamma_e;
        end
        
        function x = gamma_dot_p(obj, gamma_e, shearRate, logintMu)
        % Plastic deformation rate
            gamma_e_max = obj.gamma_e_max(logintMu);
            x = shearRate/(2 - sign(shearRate)*gamma_e/gamma_e_max);
        end
        
        function x = gamma_e_max(obj, logintMu)
        % Maximum elastic strain for a state
            d_b = 1; % backbone fractal dimension
            phi = obj.phi_a(logintMu);
            d_f = obj.par.d_f;
            x = obj.gamma_lin*((phi - obj.phi_pc)...
               ./(obj.phi_max(logintMu)- obj.phi_pc)).^((d_b+1)/(3-d_f));
        end
        
        function x = viscous_stress(obj, logintMu, shearRate)
        % Viscous stress
            x = obj.cnst.mu_s.*obj.etaTrimodal(logintMu).*shearRate;
        end
        
        function x = total_stress(obj, logintMu, gamma_e, shearRate)
        % Total stress
           x =  obj.elastic_stress(logintMu, gamma_e) ...
              + obj.viscous_stress(logintMu, shearRate);
        end
        
        function x = tau(obj, logintMu, gamma_e)
        % Relaxation time      
            x = abs(obj.gamma_lin...
                ./obj.structure_shear_rate(logintMu)...
                *abs(obj.gamma_lin/gamma_e));
        end
    end 
%% Stress responses
    methods
        out = steadyShear(obj, shearRate, initialConditions)
        
        out = steadyShearODE(obj, shearRate, initialConditions)
        
        out = stepShear(obj, initialShearRate, finalShearRate, time, init)
        
        out = UDLAOS(obj, gamma_0, omega, time, init)
        
        out = stressResponse(obj, initialShearRate, shearRate, time)
    end
%% Microstructure
    methods
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
            % Default initial guess for the 6 moments
            mn=obj.par.m_p;
            %assuming initially monodisperse
            gam=[1,1,1,1,1,1];
            K=length(gam);
            initialMu=zeros(size(gam));
            for k=1:K
            initialMu(k)=gam(k)*(mn^(k-2));
            end 
            initialMu(2)=[];
            x=log(initialMu);
            
            x = [-6.25110100853817,6.41246911055548,13.0939250626287,...
                20.1093488862100,27.4657119034575];
            
         end

         % Microstructure evolution through population balances   
         f = momicDerivative5(obj, t, shearRate,logintMu);
         
         % Elastic strain
         dX = elasticStrain(obj, ~, gamma_e, shearRate, logintMu)
         
         % Shear rate corresponding to steady shear structure 
         x = structure_shear_rate(obj, shearRate, logintMu)
   end
end