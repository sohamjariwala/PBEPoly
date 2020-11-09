load('silica.mat')
shear_rate = logspace(-1,2.3);

for i = 1:length(shear_rate)
out = steadyShear(silica, shear_rate(i));
stress(i) = out.stress;
phi_a(i) = out.phi_a;
gamma_e(i) = out.gamma_e;
end

f1 = loglog(shear_rate, stress);
    axis equal
    xlabel('Shear rate ($\dot \gamma$)',...
        'Interpreter', 'latex',...
        'FontSize', 20);
    ylabel('Shear stress ($\sigma_{xy}$)', ...
        'Interpreter', 'latex',...
        'FontSize', 20)
[LASTMSG, LASTID] = lastwarn;
if strcmp(LASTID, '')
    fprintf('\n steadyShear.m : OK PASSED\n')
else
    fprintf('\n steadyShear.m : FAILED\n')
end