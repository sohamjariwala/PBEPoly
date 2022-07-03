% clear; close all;
%% Loading the parameters and experimental data set into variables
load('silicaTest4.mat'); obj = silica;
obj = obj.inverseOfA;
%     obj.par.W = 0.008998930690468; 
%     obj.par.alfa = 0.665487794544652;
%     obj.par.b_0 = 0.003742494709682;

%     obj.par.W = 50; 
%     obj.par.alfa = 0.61;
%     obj.par.b_0 = 0.000005;
    
%shear_rate = logspace(-2, log10(300),10);
silica_SS = readmatrix('silica.ExpData/silica_SS.txt');
silica_elastic_SS = readmatrix('silica.ExpData/silica_elastic_SS.txt');
shear_rate = silica_SS(:,1);
shear_stress_SS = silica_SS(:,2);
shear_rate_elastic = silica_elastic_SS(:,1);
elastic_comp_SS = silica_elastic_SS(:,2);

stress = zeros(size(shear_rate));
logintMu = (zeros(length(shear_rate), 5));

%% Solving the steady shear equations and collecting the output variables
for i = length(shear_rate):-1:1
    if i == length(shear_rate)
        out = obj.steadyShear(shear_rate(i));
        EXITFLAG(i) = out.EXITFLAG;
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
    else
        out = obj.steadyShear(shear_rate(i), out.logintMu);
        EXITFLAG(i) = out.EXITFLAG;
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
    end
end

%% calculating Radius of gyration
c = obj.MOMIC(logintMu);
for i=1:length(shear_rate)
    phi_a(i) = obj.phi_a(logintMu(i,:));
    c(:,i) = obj.MOMIC(logintMu(i,:));
    Rg(i) =obj.fra_moment(1+1/obj.par.d_f, c(:,i))*obj.cnst.a_p*10^9*2;
    elastic_comp(i) = obj.elastic_stress(logintMu(i,:), obj.gamma_lin);
end
R_gOa_p = Rg/(obj.cnst.a_p*10^9);

error = sum(norm(stress-shear_stress_SS)./mean(shear_stress_SS))...
        /length(shear_stress_SS);
disp(error)

%% Plot generation
figure(1)
loglog(shear_rate, stress,'LineWidth',2); hold on
loglog(shear_rate, shear_stress_SS, 's','LineWidth',2);
loglog(shear_rate, elastic_comp,'LineWidth',2);
loglog(shear_rate_elastic, elastic_comp_SS, '^', 'LineWidth',2);
loglog(shear_rate, stress - elastic_comp','LineWidth',2)
xlabel('Shear rate ($\mathrm{s}^{-1}$)','Interpreter','latex','FontSize',18);
ylabel('Stress (Pa)','Interpreter','latex','FontSize',18); 
yyaxis right
plot(shear_rate, phi_a, '--','LineWidth',2);
ylabel ( '$\phi_a$','Interpreter','latex','FontSize',18);

legend('Total Shear Stress (Model)', 'Total Shear Stress (Exp.)', ...
    'Elastic Component (Model)', 'Elastic Component (Exp.)', ...
    'Viscous Component (Model)', ...
    'Volume fraction (6 moments approximation)', ...
    'Interpreter','latex','FontSize',18);

figure(2)
semilogx(shear_rate,R_gOa_p, 'LineWidth',2);
xlabel('Shear rate ($\mathrm{s}^{-1}$)','Interpreter','latex','FontSize',18);
ylabel('$R_a/a_p$', 'Interpreter','latex','FontSize',18);