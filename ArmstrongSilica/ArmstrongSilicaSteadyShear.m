%% Solving the steady shear equations and collecting the output variables
stress = zeros(size(shear_rate));
logintMu = (zeros(length(shear_rate), 5));

for i = length(shear_rate):-1:1
    if i == length(shear_rate)
        out = obj.steadyShear(shear_rate(i));
        out.A = 1;
    else
        out = obj.steadyShear(shear_rate(i), out);
    end
        EXITFLAG(i) = out.EXITFLAG;
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
        elastic_comp(i) = out.sigma_y;
        gamma_dot_p(i) = obj.gamma_dot_p(out.stress,1,out.logintMu,shear_rate(i));
        phi_max(i) = obj.phi_max(out.logintMu);
        rel_time(i) = obj.tau(out.logintMu);
end

%% Calculating radius of gyration
c = obj.MOMIC(logintMu);
for i=1:length(shear_rate)
    phi_a(i) = obj.phi_a(logintMu(i,:));
    c(:,i) = obj.MOMIC(logintMu(i,:));
    Rg(i) =obj.fra_moment(1+1/obj.par.d_f, c(:,i))*obj.cnst.a_p*10^9; %% in nm
end
R_gOa_p = Rg/(obj.cnst.a_p*10^9);

%% Plots
figure('Name', 'Steady Shear | Flow curve')
loglog(shear_rate, stress,'LineWidth',2); hold on
loglog(shear_rate, shear_stress_SS, 's','LineWidth',2);
loglog(shear_rate, elastic_comp,'LineWidth',2);
loglog(shear_rate_elastic, elastic_comp_SS, '^', 'LineWidth',2);
loglog(shear_rate, stress - elastic_comp','LineWidth',2)
xlabel('Shear rate ($\mathrm{s}^{-1}$)','Interpreter','latex','FontSize',18);
ylabel('Stress (Pa)','Interpreter','latex','FontSize',18);
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
axis([-inf inf 0.1 inf])
yyaxis right
plot(shear_rate, phi_a, '--','LineWidth',2);
ylabel ( '$\phi_a$','Interpreter','latex','FontSize',18);

legend('Total Shear Stress (Model)', 'Total Shear Stress (Exp.)', ...
    'Elastic Component (Model)', 'Elastic Component (Exp.)', ...
    'Viscous Component (Model)', ...
    'Volume fraction (6 moments approximation)', ...
    'Interpreter','latex','FontSize',10,'location','best');
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');

figure('Name', 'Steady Shear | Radius of gyration')
semilogx(shear_rate,R_gOa_p, 'LineWidth',2);
xlabel('Shear rate ($\mathrm{s}^{-1}$)','Interpreter','latex','FontSize',18);
ylabel('$R_a/a_p$', 'Interpreter','latex','FontSize',18);
axis([-inf inf -inf inf])
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
