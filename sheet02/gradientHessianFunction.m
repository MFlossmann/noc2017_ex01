
function [grad,hess]=ex2_1_function(x,y)

grad(1,1) = -2*(1-x)-400*x*(y-x*x);
grad(2,1) = 200*(y-x*x);

hess(1,1) = 2+400*(3*x*x-y);
hess(1,2) = -400*x;
hess(2,1) = -400*x;
hess(2,2) = 200;


end

