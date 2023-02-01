% 1. Steady shear flow curve
writematrix([shear_rate,stress, SSEXP.stress], ...
    "Origin_plot_data/Hippetal_Flowcurve.csv", ...
    "Delimiter",",")

% 2. Mass mean radius
writematrix([shear_rate,RG'], ...
    "Origin_plot_data/Hippetal_size_model.csv", ...
    "Delimiter",",")

writematrix([RGEXP.shear_rate, RGEXP.RG], ...
    "Origin_plot_data/Hippetal_size_experiment.csv", ...
    "Delimiter",",")