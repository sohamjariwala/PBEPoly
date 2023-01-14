%% Wei et al. Steady state figures
% 1. Flow curve
writematrix([shear_rate,stress], ...
    "Origin_plot_data/Weietal_Flowcurve_model.csv", ...
    "Delimiter",",")

% 2. Volume fraction (Structure)
writematrix([shear_rate, phi_a'], ...
    "Origin_plot_data/Weietal_Volume_fraction_model.csv", ...
    "Delimiter",",");

% 3. Transient plots
writematrix([StepDownEXP.i0p1f1_t, ...
    i1f1.stress', ...
    i2f2.stress', ...
    i3f3.stress'], ...
    "Origin_plot_data/Weietal_transient_fits_model.csv", ...
    "Delimiter",",");

writematrix([StepDownEXP.i0p1f1_t, ...
    StepDownEXP.i0p1f1_stress, ...
    StepDownEXP.i0p1f2p5_stress, ...
    StepDownEXP.i0p1f5_stress], ...
    "Origin_plot_data/Weietal_transient_fits_experiment.csv", ...
    "Delimiter",",");