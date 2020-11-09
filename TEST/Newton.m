function [x_sol iterations] = Newton( fun,jacobian,x_init,N_max,tolerance )
% Newton implements Newton's iterative method for the dolution of fun(x)=0
% INPUT:
% fun: user supplied function of x, function the root(s) of which are sought
% jacobian: user supplied function of x, the derivative of fun wrt x
% x: is the independent variable
% x_init: is the initial guess for the independent variable
% N_max: is the maximum number of iterations (10 recommended)
% tolerance: is the absolute tolerance for the solution
% OUTPUT:
% x_sol: is the solution (if converged) last guess (if not converged) 
% iterations: number of iterations (k) required for convergence within
% given tolerance;
% is equal to -(N_max+1) if no solution is found within N_max iterations;
% it takes a negative value -k if the procedure diverged after k iterations
iterations = 0;
x_sol = x_init;
diff = jacobian(x_sol)\fun(x_sol);
diff_mag_old = norm(diff);
x_sol = x_sol - diff;

for k=1:N_max
   diff = jacobian(x_sol)\fun(x_sol);
   diff_mag = norm(diff);
   if(diff_mag < tolerance)
       iterations = k;
       x_sol = x_sol - diff;
       break
   elseif(diff_mag > 10*diff_mag_old)
       iterations = -k;
       disp('Warning! Algorithm diverged!')
       break
   end
   diff_mag_old=diff_mag;
   x_sol = x_sol - diff
end
   if(iterations == 0) 
   iterations = -(N_max+1);
   disp('Warning! There is no convergence within N__max iterations')
   end
end

