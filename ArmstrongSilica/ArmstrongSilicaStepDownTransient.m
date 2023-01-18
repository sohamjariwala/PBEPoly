%% Solving transient step shear equations
initial.EXITFLAG = 1;
initial.logintMu = interp1(shear_rate, logintMu, 5);
initial.stress = interp1(shear_rate, stress,5);
initial.A = 1;
% Step down
tic; SD1 = stepShear(obj, iSD1, fSD1, time, initial); toc;
tic; SD2 = stepShear(obj, iSD2, fSD2, time, initial); toc;
tic; SD3 = stepShear(obj, iSD3, fSD3, time, initial); toc;
tic; SD4 = stepShear(obj, iSD4, fSD4, time, initial); toc;

SD1_phi_a = zeros(size(time));
SD2_phi_a = SD1_phi_a;
SD3_phi_a = SD1_phi_a;
SD4_phi_a = SD1_phi_a;

for ii = 1:length(time)
    SD1_phi_a(ii) = obj.phi_a(SD1.logintMu(ii,:));
end
for ii = 1:length(time)
    SD2_phi_a(ii) = obj.phi_a(SD2.logintMu(ii,:));
end
for ii = 1:length(time)
   SD3_phi_a(ii) = obj.phi_a(SD3.logintMu(ii,:));
end
for ii = 1:length(time)
   SD4_phi_a(ii) = obj.phi_a(SD4.logintMu(ii,:));
end

%% Plots
figure('Name', 'Step down | Stress')
box on;
semilogx(time, SD1.stress,'k-', ...
    time, SD2.stress,'r--',...
    time, SD3.stress,'m-.',...
    time, SD4.stress,'b:',...
    silica_stepDownTime, silica_stepDowni5f2p5, 'k^',...
    silica_stepDownTime, silica_stepDowni5f1p0, 'ro',...
    silica_stepDownTime, silica_stepDowni5f0p5, 'mv',...
    silica_stepDownTime, silica_stepDowni5f0p1, 'bs',...
    'MarkerSize',6,'LineWidth',2)

set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Time (s)');
ylabel(' Stress (Pa)');
legend('2.5 s^{-1}','1 s^{-1}','0.5 s^{-1}','0.1 s^{-1}');
axis([-inf inf 0 25])
grid on;

figure('Name', 'Step down | Volume fraction')
box on;
semilogx(time, SD1_phi_a,'k-', ...
        time, SD2_phi_a,'r--',...
        time, SD3_phi_a,'m-.',...
        time, SD4_phi_a,'b:',...
        'MarkerSize',6,'LineWidth',2)

    set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Time (s)');
ylabel('\phi_a');
legend('2.5 s^{-1}','1 s^{-1}','0.5 s^{-1}', '0.1 s^{-1}');
axis([0.001 tEnd -inf inf])
grid on;