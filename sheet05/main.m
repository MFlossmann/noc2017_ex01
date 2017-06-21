clear all;

import casadi.*

N = 50;
nx = 1;
nu = 1;
x0 = 0.5;
h = 0.1;

x = MX(N+1,nx);
u = MX.sym('u',N,nu);
n = MX.sym('n',1);
Uk = ones(N,nu);

%for i=1:N-1
    %u(i) = 1;
%end

% Specify phi
x(1) = x0;
for i=1:N
    x(i+1) = x(i) + h*( (1-x(i) ) * x(i) + u(i)); 
end

Phi = Function('Phi', {u}, {x});
J = Function('J', {u}, {jacobian(x,u)});
 
%% b) Forward AD

AD_diff = forw_AD(Uk, x0, nx, h);
CA_diff = full(J(Uk));

[CA_diff(end,:)', AD_diff(end,:)']

%% c) Backward AD


%% d) J_FAD

clear AD_diff;

m= 5;

AD_diff = J_FAD(Uk, m, x0, nx, h);
CA_diff = full(J(Uk));


tmp = size(CA_diff,1) - m;
CA_diff = CA_diff(tmp:end,:);

errors = ones(m,1);
for i=1:m
    errors(i) = norm(AD_diff(i,:) - CA_diff(i,:))
end

display(errors)