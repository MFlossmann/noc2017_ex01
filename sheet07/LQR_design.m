function LQR_design(Q,R,Ts,nSteps)
% linearize
input.Ts = Ts; input.nSteps = nSteps;
x_lin = [0; 0]; u_lin = 0;
input.x = x_lin; input.u = u_lin;
output = RK4_integrator( @ode, input );
A = output.sensX;
B = output.sensU;

% design LQR controller
% INSERT YOUR CODE HERE:
%[K,P] = ...;
P = Q;
tmp = (R + B'*P*B);
K = tmp\(B'*P*A);

[K, S] = dlqr(A, B, Q, R, 0);


save lqr.mat A B Q R K P
end