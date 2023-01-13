%% Transient step shear plot
initial.EXITFLAG = 1;
initial.logintMu = interp1(SSEXP.shear_rate, logintMu, 0.1);
initial.stress = interp1(SSEXP.shear_rate, stress, 0.1);
initial.A = 1;
%% Set Step Shear parameters
i1 = 0.1; f1 = 1;
i2 = 0.1; f2 = 2.5;
i3 = 0.1; f3 = 5;

%% Solution
i1f1 = stepShear(obj, i1, f1, StepDownEXP.i0p1f1_t, initial);
i2f2 = stepShear(obj, i2, f2, StepDownEXP.i0p1f2p5_t, initial);
i3f3 = stepShear(obj, i3, f3, StepDownEXP.i0p1f5_t, initial);


%% Error
transient_error = (norm((i1f1.stress-StepDownEXP.i0p1f1_stress)./mean(StepDownEXP.i0p1f1_stress)) + ...
                   norm((i2f2.stress-StepDownEXP.i0p1f2p5_stress)./mean(StepDownEXP.i0p1f2p5_stress)) + ...
                   norm((i3f3.stress-StepDownEXP.i0p1f5_stress)./mean(StepDownEXP.i0p1f5_stress)))...
                 ./length(StepDownEXP.i0p1f5_stress)/3;

fprintf("Transient step shear error = %f\n", transient_error);

 %% Transient plots
 figure('Name','Step down transients')
    semilogx(StepDownEXP.i0p1f1_t, StepDownEXP.i0p1f1_stress,'ko',...     
             StepDownEXP.i0p1f2p5_t, StepDownEXP.i0p1f2p5_stress, 'bv', ...
             StepDownEXP.i0p1f5_t, StepDownEXP.i0p1f5_stress, 'g^',...
             StepDownEXP.i0p1f5_t, i3f3.stress, 'g',...
             StepDownEXP.i0p1f1_t, i1f1.stress,'k',...
             StepDownEXP.i0p1f2p5_t, i2f2.stress, 'b', ...
            'MarkerSize',6,'LineWidth',2)
         
             legend('1 (1/s)','2.5 (1/s)','5 (1/s)');
    set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
    xlabel('Time (s)');
    ylabel(' Stress (Pa)');

    axis([-inf inf 0 12])
    grid on;