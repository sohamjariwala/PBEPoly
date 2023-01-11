obj = PBEPoly;

%% Steady state
shear_rate = logspace(-2,2,20);

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

loglog(shear_rate, stress, ...
    shear_rate, elastic_stress, ...
    'LineWidth', 2);
xlabel('$\dot\gamma$','Interpreter','latex');
ylabel('$\sigma_{yx}s$','Interpreter','latex');
set(gca,'FontSize',20)