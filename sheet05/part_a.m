
clear all;

import casadi.*;

x_0 = 0.5;
N = 50;
h = 0.1;
u_val = SX.sym('u_val',N,1);
u_val(:) = 1;
u = SX.sym('u',N,1);
x = SX.sym('x',N,1);
x(1) = x_0;

for i = 2 : N
    x(i) = x(i-1)+h*((1-x(i-1))*x(i-1)+u(i));
end

Phi = Function('Phi',{u},{x});
JPhi = jacobian(x,u);
jacobian_Phi = Function('jacobian_Phi',{u},{JPhi});

display(Phi(u_val));
display(jacobian_Phi(u_val));
