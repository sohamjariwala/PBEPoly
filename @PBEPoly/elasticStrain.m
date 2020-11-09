function dX = elasticStrain(obj, ~, gamma_e, shearRate, logintMu)
    %gamma_e = min([gamma_e, obj.gamma_lin]);
    dX = (shearRate-gamma_e/obj.tau( logintMu, gamma_e));
end