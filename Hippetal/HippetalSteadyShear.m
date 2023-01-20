% shear_rate = logspace(-2,2,20)';
shear_rate = SSEXP.shear_rate;
stress = zeros(size(shear_rate));
logintMu = (zeros(length(shear_rate), 5));

%% Solving the steady shear equations and collecting the output
for i = length(shear_rate):-1:1
    if i == length(shear_rate)
        out = obj.steadyShear(shear_rate(i));
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
        phi_a(i) = out.phi_a;
        phi_max(i) = obj.phi_max(out.logintMu);
    else
        out = obj.steadyShear(shear_rate(i), out);
        out = obj.steadyShearODE(shear_rate(i), out);
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
        phi_a(i) = out.phi_a;
        phi_max(i) = obj.phi_max(out.logintMu);
    end
end

%% Calculating radius of gyration
c = obj.MOMIC(logintMu);
for i=1:length(shear_rate)
    phi_a(i) = obj.phi_a(logintMu(i,:));
    c(:,i) = obj.MOMIC(logintMu(i,:));
    RG(i) =obj.fra_moment(1+1/obj.par.d_f, c(:,i))*obj.cnst.a_p; %% in m
end

%% Steady shear plot
% Plot flow curve and volume fraction
figure('Name','Steady state flow curve')
loglog(shear_rate, stress, SSEXP.shear_rate, SSEXP.stress,'o',...
    'MarkerSize',6,'LineWidth',2)
xlabel('Shear rate (s^{-1})','FontSize',18);
ylabel('Stress (Pa)','FontSize',18);
axis([-inf inf -inf inf]);
yyaxis right
loglog(shear_rate,phi_a,'LineWidth',2);
ylabel('Volume fraction (-)')
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');

% Plot radius of gyration
figure('Name', 'Steady Shear | Radius of gyration')
semilogx(shear_rate, RG/min(RG),RGEXP.shear_rate, RGEXP.RG'/min(RGEXP.RG), 's', 'LineWidth',2);
xlabel('Shear rate ($\mathrm{s}^{-1}$)','Interpreter','latex','FontSize',18);
ylabel('$R_a$ (m)' , 'Interpreter','latex','FontSize',18);
axis([-inf inf -inf inf])
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');

% Plot distribution
obj.distribution(shear_rate,logintMu);