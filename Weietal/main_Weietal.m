addpath('../');
%-------------------------------------------------------------------------------
parVec = [3.091690100862115   0.604932033000000  -6.177390351405440   2.300037015000000   0.879680247000000   5.732836850990885  36.539401529999999   0.795985449000000   0.323856153000000];
obj = PBEPoly;
%% Steady state
    %% Loading the parameters 
    obj.par.W = exp(parVec(1))-1;
    obj.par.alfa = parVec(2);
    obj.par.b_0 = exp(parVec(3));
    obj.par.d_f = parVec(4);
    obj.par.porosity = parVec(5);
    obj.par.m_p = exp(parVec(6));
    obj.par.p = 4;

    obj.cnst.G_0 = parVec(7);
    obj.cnst.sigma_y0 = parVec(8);
    obj.cnst.mu_s = parVec(9);
    

loadExperimentalData;

WeietalSteadyShear;

WeietalStepDownTransient;

obj.distribution(shear_rate,logintMu);

writeDataToOrigin

%% Error
fObj = SS_error + transient_error;
fprintf("Total error = %f\n", fObj);