shear_rate = SSEXP.shear_rate;
stress = zeros(size(SSEXP.shear_rate));
logintMu = (zeros(length(SSEXP.shear_rate), 5));

%% Solving the steady shear equations and collecting the output
for i = length(SSEXP.shear_rate):-1:1
    if i == length(SSEXP.shear_rate)
        out = obj.steadyShear(SSEXP.shear_rate(i));
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
        phi_a(i) = out.phi_a;
    else
        out = obj.steadyShear(SSEXP.shear_rate(i), out);
        out = obj.steadyShearODE(SSEXP.shear_rate(i), out);
        stress(i) = out.stress;
        logintMu(i,:) = out.logintMu;
        phi_a(i) = out.phi_a;
    end
end

SS_error = norm((stress-SSEXP.stress)./(SSEXP.stress))...
    /length(SSEXP.stress);

fprintf("Steady state error = %f\n", SS_error);

%% Steady shear plot
figure('Name','Steady state flow curve')
loglog(SSEXP.shear_rate, stress, SSEXP.shear_rate, SSEXP.stress,'o',...
    'MarkerSize',6,'LineWidth',2)
xlabel('Shear rate (s^{-1})','FontSize',18);
ylabel('Stress (Pa)','FontSize',18);
axis([-inf inf 5e-1 100]);
yyaxis right
loglog(SSEXP.shear_rate,phi_a,'LineWidth',2);
ylabel('Volume fraction (-)')
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
