function dX = elasticStrain(obj, ~, gamma_e, shearRate, logintMu)
gamma_e = max([0, min([gamma_e, obj.gamma_lin])]);
gamma_max = obj.gamma_lin;
    
dX = shearRate*(1 - obj.structure_shear_rate(shearRate,logintMu)/shearRate*(gamma_e/gamma_max)^4);

end