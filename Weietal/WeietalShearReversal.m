initial.EXITFLAG = 1;
initial.logintMu = interp1(shear_rate, logintMu, 0.75);
initial.stress = -interp1(shear_rate, stress, 0.75);
initial.A = -1;

time = logspace(-5,log10(5000));
rev1 = obj.shearReversal(-0.75, time, initial);

figure('Name', 'Shear flow reversal')
semilogx(time, rev1.stress, 'LineWidth',2)
xlabel('Time (s)')
ylabel("\sigma/\sigma_\infty")
axis([-inf inf -inf inf])
set(gca, 'FontName', 'Times', 'FontSize', 20, 'LineWidth', 2)