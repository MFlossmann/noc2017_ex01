function [Tout, Xout] = rk4_CasADi(ode, h, numSteps, x0)

    Tout = 0:h:h*(numSteps-1);
    
% Declare model variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1; x2];

% Model equations
xdot = ode(0,X);
f = Function('IVP', {X}, {xdot});

% Formulate RK4 integrator
X0 = MX.sym('X0', 2);
X = X0;

k1 = f(X);
k2 = f(X + h/2 * k1);
k3 = f(X + h/2 * k2);
k4 = f(X + h * k3);
X = X + h/6*(k1 + 2*k2 + 2*k3 + k4);
F = Function('F', {X0}, {X}, {'x0'}, {'x'});

    % Run simulation
     Xout = zeros(length(X0),numSteps);
     Xout(:,1) = x0;
     for k = 2:1:numSteps
        Fk = F('x0',Xout(:,k-1));
        Xout(:,k) = full(Fk.x);
     end
end