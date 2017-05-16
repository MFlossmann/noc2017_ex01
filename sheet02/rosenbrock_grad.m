function [ grad ] = rosenbrock_grad( x, y )
%rosenbrock_grad Returns the gradient of the Rosenbrock problem
    grad = zeros(2,1);

    grad(1) = 2*(x*(200*(x^2-y) + 1) - 1);
    grad(2) = 200*(y-x^2);
end

