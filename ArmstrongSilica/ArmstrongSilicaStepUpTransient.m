% Step up
initial.EXITFLAG = 1;
initial.logintMu = interp1(shear_rate, logintMu, 0.1);
initial.stress = interp1(shear_rate, stress,0.1);

% Step up
tic; SU1 = stepShear(obj, iSU1, fSU1, time, initial); toc;
tic; SU2 = stepShear(obj, iSU2, fSU2, time, initial); toc;
tic; SU3 = stepShear(obj, iSU3, fSU3, time, initial); toc;
tic; SU4 = stepShear(obj, iSU4, fSU4, time, initial); toc;
tic; SU5 = stepShear(obj, iSU5, fSU5, time, initial); toc;

SU1_phi_a = zeros(size(time));
SU2_phi_a = SU1_phi_a;
SU3_phi_a = SU1_phi_a;
SU4_phi_a = SU1_phi_a;
SU5_phi_a = SU1_phi_a;

for ii = 1:length(time)
   SU1_phi_a(ii) = obj.phi_a(SU1.logintMu(ii,:));
end
for ii = 1:length(time)
   SU2_phi_a(ii) = obj.phi_a(SU2.logintMu(ii,:));
end
for ii = 1:length(time)
   SU3_phi_a(ii) = obj.phi_a(SU3.logintMu(ii,:));
end
for ii = 1:length(time)
   SU4_phi_a(ii) = obj.phi_a(SU4.logintMu(ii,:));
end
for ii = 1:length(time)
   SU5_phi_a(ii) = obj.phi_a(SU5.logintMu(ii,:));
end

%% Plots
figure('Name', 'Step up | Stress')
box on;
semilogx(time, SU1.stress,'k-', ...
    time, SU2.stress,'r--',...
    time, SU3.stress,'m-.',...
    time, SU4.stress,'b:',...
    time, SU5.stress,'g.-',...
    silica_stepUpTime, silica_stepUpi0p1f5, 'k^',...
    silica_stepUpTime, silica_stepUpi0p1f2p5, 'ro',...
    silica_stepUpTime, silica_stepUpi0p1f1, 'mv',...
    silica_stepUpTime, silica_stepUpi0p1f0p5, 'bs',...
    silica_stepUpTime, silica_stepUpi0p1f0p25, 'gp',...
    'MarkerSize',6,'LineWidth',2)

set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Time (s)');
ylabel(' Stress (Pa)');
legend('5 s^{-1}','2.5 s^{-1}','0.1 s^{-1}','0.5 s^{-1}','2.5 s^{-1}');
axis([-inf inf 0 inf])
grid on;

figure('Name', 'Step up | Volume fraction')
box on;
semilogx(time, SU1_phi_a,'k-', ...
        time, SU2_phi_a,'r--',...
        time, SU3_phi_a,'m-.',...
        time, SU4_phi_a,'b:',...
        time, SU5_phi_a,'g.-',...
        'MarkerSize',6,'LineWidth',2)

    set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
xlabel('Time (s)');
ylabel('\phi_a');
legend('5 s^{-1}','2.5 s^{-1}','1.0 s^{-1}','0.5 s^{-1}','0.25 s^{-1}');
axis([-inf inf -inf inf])
grid on;