function [ hess ] = rosenbrock_hessian( x,y )
%rosenbrock_hessian Calculates the hessian of the Rosenbrock problem
    hess = zeros(2,2);
    
    hess(1,1) = 800*x + 2;
    hess(1,2) = -400;
    hess(2,1) = hess(1,2);
    hess(2,2) = 200;
end

