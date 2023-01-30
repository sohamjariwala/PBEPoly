rev_count = 0;
for rate = [0.1 0.5 2.5 5 10 20]
    initial.EXITFLAG = 1;
    initial.logintMu = interp1(shear_rate, logintMu, rate);
    initial.stress = -interp1(shear_rate, stress, rate);
    initial.A = -1;

    strain = logspace(-5,log10(100));
    rev_count = rev_count + 1;
    rev(rev_count) = obj.shearReversal(-rate, strain/rate, initial);
end