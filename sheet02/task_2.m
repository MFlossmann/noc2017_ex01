% Numerical Optimal control Ex. Sheet 2, Task 2

clear variables
%close all
import casadi.*

v = SX.sym('v',4,1);

% Define the position through the function
tic
p_T = ballistic_dynamics_RK4(v);
toc

% J is the distance between the projectiles
J = norm(p_T(1:2) - p_T(3:4));
%%
% Restrictions: the z coordinates have to be > 0
g = [norm(v(1:2)),... % maximum velocity amount
     norm(v(3:4)),... % 
     p_T(2),...
     p_T(4)];
%g = [];

nlp = struct('x', v, 'f', J', 'g', g);

solver = nlpsol('solver', 'ipopt', nlp);

res = solver('x0' , [3, 10, 3, -10],...       % solution guess
             'lbg', [0,0,0,0],...           % lower bound on g
             'ubg', [15,15,inf,inf],...     % upper bound on g
             'lbx', [-inf; 0; -inf; 0],...  % lower bound on x
             'ubx', [inf; inf; inf; inf]);  % upper bound on x
         
v_opt = full(res.x);
fprintf('Optimal speed:');
display(v_opt);
p_opt = ballistic_dynamics_RK4(v_opt);
distance_opt = norm(p_opt(1:2) - p_opt(3:4));

fprintf('End distance: %s\n',num2str(distance_opt));

%%
M = 100;
T = 0.5;
DT = T/M;

% simulation
X0 = [0, 0, v_opt(1), v_opt(2), 10, 0, v_opt(3), v_opt(4)];
X_sim = zeros(8,M);
X_sim(:,1) = X0;
for i=2:M
    k_1 = ode(X_sim(:,i-1));
    k_2 = ode(X_sim(:,i-1) + k_1*DT/2);
    k_3 = ode(X_sim(:,i-1) + k_2*DT/2);
    k_4 = ode(X_sim(:,i-1) + k_3*DT);
    
    X_sim(:,i) = X_sim(:,i-1) + DT * (k_1 + 2*k_2 + 2*k_3 + k_4) / 6;
end

figure;
hold on;
plot(X_sim(1,:),X_sim(2,:));
plot(X_sim(5,:),X_sim(6,:));
hold off;

%% BONUS
% no drag
p_T = ballistic_dynamics_RK4(v, 0, 0);

% J is the distance between the projectiles
J = norm(p_T(1:2) - p_T(3:4));

% Restrictions: the z coordinates have to be > 0
g = [p_T(2); p_T(4)];

nlp = struct('x', v, 'f', J', 'g', g);

solver = nlpsol('solver', 'ipopt', nlp);

res = solver('x0' , [5, 5, -5, 5],...   % solution guess
             'lbg', [0;0],...           % lower bound on g
             'lbx', -inf,...            % lower bound on x
             'ubx',  inf);              % upper bound on x
         
v_opt = full(res.x);
fprintf('Optimal speed:');
display(v_opt);
p_opt = ballistic_dynamics_RK4(v_opt, d1, d2);
distance_opt = norm(p_opt(1:2) - p_opt(3:4));

fprintf('End distance: %s\n',num2str(distance_opt));

%% Simulation

M = 100;
T = 0.5;
DT = T/M;

d1 = 0;
d2 = 0;

% simulation
X0 = [0, 0, v_opt(1), v_opt(2), 10, 0, v_opt(3), v_opt(4)];
X_sim = zeros(8,M);
X_sim(:,1) = X0;
for i=2:M
    k_1 = ode(X_sim(:,i-1), d1, d2);
    k_2 = ode(X_sim(:,i-1) + k_1*DT/2, d1, d2);
    k_3 = ode(X_sim(:,i-1) + k_2*DT/2, d1, d2);
    k_4 = ode(X_sim(:,i-1) + k_3*DT, d1, d2);
    
    X_sim(:,i) = X_sim(:,i-1) + DT * (k_1 + 2*k_2 + 2*k_3 + k_4) / 6;
end

%%
figure;
hold on;
plot(X_sim(1,:),X_sim(2,:));
plot(X_sim(5,:),X_sim(6,:));
hold off;