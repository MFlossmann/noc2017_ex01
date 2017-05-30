% Implementation of an interior point method 
clear variables
close all
clc

import casadi.*
% Problem definition
nv = 2;
ne = 1;
ni = 1;

DELTA = 1e-3;
x_test = [2,3];

x = MX.sym('x',nv);

f = (x(1)- 4)^2 + (x(2) - 4)^2;
g = sin(x(1)) - x(2)^2;
h = x(1)^2 + x(2)^2 - 4;


% YOUR CODE HERE:
% Inequality contstraints
H = Function('H',{x},{h});
h_test = 9.0;
if full(H(x_test)) - h_test > DELTA
    error('Failure! H doesn`t fit!');
end

% Equality contstraints
G = Function('G',{x},{g});
G_test = -8.0907;
if full(G(x_test)) - G_test > DELTA
    error('Failure! G doesn`t fit!');
end
%
% Objective
F = Function('F',{x},{f});
F_test = 5;
if full(F(x_test)) - F_test > DELTA
    error('Failure! F doesn`t fit!');
end

% Jacobian of g
Jg = Function('Jg',{x},{jacobian(g,x)});
Jg_test = [-0.4161,   -6.0000];
if full(Jg(x_test)) - Jg_test > DELTA
    error('Failure! Jg doesn`t fit!');
end

% Jacobian of h
Jh = Function('Jh',{x},{jacobian(h,x)});
Jh_test = [4,   6];
if full(Jh(x_test)) - Jh_test > DELTA
    error('Failure! Jh doesn`t fit!');
end

% Jacobian of f
Jf = Function('Jf',{x},{jacobian(f,x)});
Jf_test = [-4,   -2];
if full(Jf(x_test)) - Jf_test > DELTA
    error('Failure! Jf doesn`t fit!');
end

% Hessian of g
Hg = Function('Hg',{x},{hessian(g,x)});
Hg_test =    [  -0.9093,  0        ; ...
                0,        -2.0000 ];
if full(Hg(x_test)) - Hg_test > DELTA
    error('Failure! Hg doesn`t fit!');
end            
% Hessian of h
Hh = Function('Hh',{x},{hessian(h,x)});
Hh_test =   [   2, 0; ...
                0  2];
if full(Hh(x_test)) - Hh_test > DELTA
    error('Failure! Hh doesn`t fit!');
end

% Hessian of f
% Hf = Function(...);
Hf = Function('Hf',{x},{hessian(f,x)});
Hf_test =   [   2, 0; ...
                0  2];
if full(Hf(x_test)) - Hf_test > DELTA
    error('Failure! Hf doesn`t fit!');
end
            
% Interior point solver
max_it = 100;
xk = [-2;-4];
lk = 10*ones(ne,1);
vk = 10*ones(ni,1);
sk = 10*ones(ni,1);
iter = zeros(nv + ne + ni + ni,max_it);
iter(:,1) = [xk;lk;vk;sk];
tau = 2;
k_b = 1/3;
th_1 = 1.0e-8;
th_2 = 1.0e-8;
for i = 2:max_it
    % Build KKT system
    Hf_e    = Hf(xk);
    Hg_e    = Hg(xk);
    Hh_e    = Hh(xk);
    Hl      = full(Hf_e) + full(Hg_e)*lk + full(Hh_e)*vk;
    
    Jg_e    = Jg(xk);
    Jh_e    = Jh(xk);
    Jf_e    = Jf(xk);
    
    g_e     = G(xk);
    h_e     = H(xk);
    
    % YOUR CODE HERE:
    % Buiild the KKT system
    M = full([Hf_e + Hg_e*lk + Hh_e*vk, Jg_e', Jh_e', zeros(2,1);...
              Jg_e, 0, 0, 0; ...
              Jh_e, 0, 0, 1; ...
              0, 0, 0, sk, vk]);  
    
    rhs = - full([    (Jf_e + Jg_e*lk + Jh_e*vk)'; ...
                      g_e; ...
                      h_e + sk; ...
                      sk*vk - tau]);
   
    % Termination condition
    if norm(rhs) < th_1
        if tau < th_2
            display('Solution found!')
            break;
        else
            tau = tau*k_b;
        end
    end
     
    % YOUR CODE HERE:
    % Compute Newton step
    sol = inv(M)*rhs;
    
    % line-serach
    max_ls = 100;
    x_step  = sol(1:nv);
    l_step  = sol(nv+1:nv+ne);
    v_step  = sol(nv+ne+1:nv+ne+ni);
    s_step  = sol(nv+ne+ni+1:end);
    alpha = 1;
    k_ls = 0.9;
    min_step = 1.0e-8;
    for j=1:max_ls
        
        % YOUR CODE HERE: 
        % Compute trial step
        v_t = vk + alpha*v_step;
        s_t = sk + alpha*s_step;
        if (isempty(v_t(v_t <= 0)) && isempty(s_t(s_t <= 0)))
            break;
        end
        
        % YOUR CODE HERE:
        % Decrease alpha
        
        % 0.7 seems like a good value. 0.9 doesn't converge and 0.5
        % terminates
        alpha = alpha*0.9;
        
        % Terminiation condition
        if norm(alpha*[ v_step;s_step]) < min_step
            error('Line search failed! Could not find dual feasible step.')
        end
    end
    
    xk  = xk + alpha*x_step;
    lk  = lk + alpha*l_step;
    vk  = vk + alpha*v_step;
    sk  = sk + alpha*s_step;
    
    
    % Print some results
    display(['Iteration: ',num2str(i)]);
    display(['tau = ',num2str(tau)]);
    display(['norm(rhs) = ',num2str(norm(rhs))]);
    display(['step size = ',num2str(alpha)])
    iter(:,i) = [xk;lk;vk;sk];
end
iter = iter(:,1:i-1);
plot(iter')
grid on
xlabel('iterations')
ylabel('solution')

% Plot feasible set, and iterations
figure()
pl = ezplot('sin(x) - y^2');
set(pl,'Color','red');
hold all
pl= ezplot('x^2 + y^2 - 4');
set(pl,'Color','blue');
ezcontour('(x- 4)^2 + (y - 4)^2')
plot(iter(1,:),iter(2,:),'--','Color','black')
plot(iter(1,:),iter(2,:),'o','Color','black')
title('Iterations in the primal space')
grid on