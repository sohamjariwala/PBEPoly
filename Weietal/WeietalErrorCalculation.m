%% Steady shear error 
SS_error = norm((stress-SSEXP.stress)./(SSEXP.stress))...
    /length(SSEXP.stress);

fprintf("Steady state error = %f\n", SS_error);

%% Transient Error
transient_error = (norm((i1f1.stress-StepDownEXP.i0p1f1_stress)./mean(StepDownEXP.i0p1f1_stress)) + ...
                   norm((i2f2.stress-StepDownEXP.i0p1f2p5_stress)./mean(StepDownEXP.i0p1f2p5_stress)) + ...
                   norm((i3f3.stress-StepDownEXP.i0p1f5_stress)./mean(StepDownEXP.i0p1f5_stress)))...
                 ./length(StepDownEXP.i0p1f5_stress)/3;

fprintf("Transient step shear error = %f\n", transient_error);


transient_error_test = norm(i4f4.stress - StepDownEXP_Test.i3f9_stress)./mean(StepDownEXP_Test.i3f9_stress) ...
    /length(StepDownEXP_Test.i3f9_stress);


%% Total error
fObj = SS_error + transient_error;
fprintf("Total error = %f\n", fObj);