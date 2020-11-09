classdef PBEPoly
% A rheological constitutive model for carbon black via population
% balances using monodisperse closure.
    
%% Parameters    
    properties 
       % Parameters to be determined from fitting  
       par = struct(...
            'W', [], ...                % Fuch's stability ration
            'alfa', [], ...             % Collision efficiency
            'b_0' , [], ...             % Breakage constant
            'd_f' , [], ...             % Fractal dimension
            'porosity', [], ...         % R_h/R_a
            'm_p', [], ...              % #(particles) in a primary cluster
            'bexp', [] ...
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
            x = obj.cnst.phi_p*exp(logintMu(1))^(1-3/obj.par.d_f);
        end
                
        function x = cutOff(obj, logintMu)
        % Hyperbolic cutoff function to model dynamic arrest
            phi_a = obj.phi_a(logintMu);
            x = tanh(2.65 * (obj.phi_max(logintMu) - phi_a)...
                        ./ (obj.phi_max(logintMu) - obj.phi_pc));
        end
        
        function x = phi_max(obj, logintMu)
        % Maximum packing fraction as a function of moments (polydisperse)
            c=obj.MOMIC(logintMu);
            phi_rcp = obj.cnst.phi_max;
            df = obj.par.d_f;
            
            fra_moment = @(j,c) obj.fra_moment(j,c);
            
            x = phi_rcp* ...
                exp(0.2* ...
                sqrt( log( ...
                fra_moment(4/df,c)*fra_moment(2/df,c)/fra_moment(3/df,c)...
                )));
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
        x = obj.cnst.G_0.*((obj.phi_a(logintMu) - obj.cnst.phi_p)...
         ./(obj.phi_max(logintMu) - obj.cnst.phi_p)).^(2/(3-obj.par.d_f));
        end
        
        function x = elastic_stress(obj, logintMu, gamma_e)
        % Elastic stress
            x = obj.G(logintMu).*gamma_e;
        end
        
        function x = viscous_stress(obj, logintMu, shear_rate)
        % Viscous stress
            x = obj.cnst.mu_s.*obj.etaTrimodal(logintMu).*shear_rate;
        end
        
        function x = total_stress(obj, logintMu, gamma_e, shear_rate)
        % Total stress
           x =  obj.elastic_stress(logintMu, gamma_e) ...
              + obj.viscous_stress(logintMu, shear_rate);
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
        out = steadyShear(obj, shear_rate)
        
        out = stepShear(obj, initialShearRate, finalShearRate, time)
        
        out = UDLAOS(obj, gamma_0, omega, time)
        
        out = stressResponse(obj, initialShearRate, shearRate, time)
    end
%% Microstructure
    methods
        % MOMIC
        function x = MOMIC(obj, logintMu)
        %Using MOMIC Method to find the coefficients from integer moments
            x=obj.lnvA*(logintMu');
        end
        
         function x = get.lnvA(~)
         % Create the A matrix for 
            K=9;
            A=zeros(K,K);
            j=[0,2,3,4,5,6,7,8,9];
            k=j;
                for a=1:K
                    for b=1:K
                         A(a,b)=j(a)^k(b)-j(a);
                    end 
                end  
            x=inv(A);
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
         
         function x=InitialCondition(~)
            % Default initial guess for the 9 moments
            mn=600;
            %assuming initially monodisperse
            gam=[1,1,1,1,1,1,1,1,1,1];
            K=length(gam);
            initialMu=zeros(size(gam));
            for k=1:K
            initialMu(k)=gam(k)*(mn^(k-2));
            end 
            initialMu(2)=[];
            x=log(initialMu);
         end

         % Microstructure evolution through population balances   
         f = momicDerivative9(obj, t, shearRate,logintMu);
         
         % Elastic strain
         dX = elasticStrain(obj, ~, gamma_e, shearRate, logintMu)
         
         % Shear rate corresponding to steady shear structure 
         x = structure_shear_rate(obj, phi_a, shear_rate)
   end
end