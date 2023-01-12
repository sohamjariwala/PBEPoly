%% Distribution generation function for MOMIC
% This function generates lognormal distribution for provided physical
% moments in log scale.
% INPUT: Log scale moments
% OUTPUT: distribution function and plot

time = shear_rate;
c=obj.MOMIC(logintMu);

gamma2 = exp(logintMu(:,1)).*exp(logintMu(:,2));

a_p = 8e-9;
d_f = 2.11;
N = 1000;

figure;c_map = colormap(winter(length(time)));
lnsigmag = zeros(size(time));
mulnX = zeros(size(time));
for i = 1:length(time)
  lnsigmag(i) = sqrt(log(gamma2(i)));
  mulnX(i) = 1/sqrt(gamma2(i));
  average(i) = exp(mulnX(i) + lnsigmag(i)^2/2);
  a = linspace(-8,8,N);
  b =  log(10)*(exp(a).^(1/d_f))...
      .*(1/(sqrt(2*pi)*lnsigmag(i))...
      *exp(-0.5*((a - log(mulnX(i)))/(lnsigmag(i))).^2));
  a = 10^6*(exp(a)*obj.fra_moment(1+1/d_f,c(:,i))*2*a_p);
  semilogx(a', b', 'Color', c_map(i,:), 'LineWidth',2);
  hold on
end

xlabel('Aggregate size (diameter) [\mum]');
ylabel('ln(10) x^{1+1/d_f} f(x)');
title('Probability density vs. size')
axis([-inf, inf, 0, inf])
axis square
clim([shear_rate(1) shear_rate(end)])
set(gca,'FontSize',14,'FontWeight','bold','linewidth',2, 'FontName','Times');
colorbar, grid on