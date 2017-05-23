clear all
close all
clc

import casadi.*

PLOT_GRADIENTS = 1;

T = 0.5; % time horizon
M = 100; % discretization steps

% Declare model variables
x = MX.sym('x',4,1);

alpha_min =  -1;

if PLOT_GRADIENTS
    n_points = 100;
    alpha_max = 1;
else
    n_points = 1;
    alpha_max = -alpha_min;
end

alpha_v = linspace(alpha_min,alpha_max,n_points);
alpha = MX.sym('alpha',1,1);

% Initial position of the balls
x01 = 0;
x02 = 10;

% max velocity
v_max = 18;

% Model equations
xdot = ballistic_dynamics(x);

% Continuous time dynamics
f = Function('f', {x}, {xdot});

%% Formulate discrete time dynamics
% Fixed step Runge-Kutta 4 integrator
DT = T/M;
X0 = MX.sym('X0', 2, 1);

X = [x01; 0; X0];
for j=1:M
    k1 = f(X);
    k2 = f(X + DT/2 * k1);
    k3 = f(X + DT/2 * k2);
    k4 = f(X + DT * k3);
    X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
end

F = Function('F', {X0}, {X});

% A second function for plotting is defined below 
X0_plot = MX.sym('X0_plot', 4, 1);
DTs = MX.sym('DTs', 1, 1);

X = X0_plot;
M_plot = 10;
for j=1:M_plot
    k1 = f(X);
    k2 = f(X + DTs/2 * k1);
    k3 = f(X + DTs/2 * k2);
    k4 = f(X + DTs * k3);
    X=X+DTs/6*(k1 +2*k2 +2*k3 +k4);
end

Fs = Function('Fs', {X0_plot,DTs}, {X});

%% NLP Formulation
% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% "Lift" initial conditions
V01 = MX.sym('V01', 1, 1);
V02 = MX.sym('V02', 1, 1);

% Insert you code here (NLP formulation)
% ...
% ...

% Insert you code here (Jacobian of the constraints)
% ...
% ...
% ...

% Insert you code here (create an NLP solver)
% ...
% ...
% ...

gradients = zeros(6,n_points);
lambdas = zeros(3,n_points);

h = figure(1);
subplot(211)

hold all
p1 = plot([0, 0], [0, 0]);
p2 = plot([0, 0], [0, 0]);
p3 = plot([0, 0], [0, 0]);
p4 = plot(0,0);

for i = 1:n_points
    
    alpha_val = alpha_v(i);
    
    % Insert you code here (solve the NLP)
    % ...
    % ...
    
    w_opt = full(sol.x);
    Jg_eval = full(Jg(w_opt,alpha_val));
    gradients(:,i) = reshape(Jg_eval.',6,1); 
    
    if PLOT_GRADIENTS
        s1 = subplot(211);
        x0=10;
        y0=10;
        width=400;
        height=500;
        set(gcf,'units','points','position',[x0,y0,width,height])
        
        delete(p1)
        delete(p2)
        delete(p3)
        
        p1 = plot([0, Jg_eval(1,1).'/norm(Jg_eval(1,:))], [0, Jg_eval(1,2).'/norm(Jg_eval(1,:))],'r','Linewidth',2);
        p2 = plot([0, Jg_eval(2,1).'/norm(Jg_eval(2,:))], [0, Jg_eval(2,2).'/norm(Jg_eval(2,:))],'b','Linewidth',2);
        p3 = plot([0, Jg_eval(3,1).'/norm(Jg_eval(3,:))], [0, Jg_eval(3,2).'/norm(Jg_eval(3,:))],'g','Linewidth',2);
        
        xlabel('y')
        ylabel('z')
        grid on
        
        xlim([-1,1])
        ylim([-1,1])
        
        lambdas(:,i) = full(sol.lam_g);
        
        subplot(212)
        delete(p4)
        
        p4 = plot(lambdas(:,1:i).','Linewidth',1);
        xlabel('iterations')
        ylabel('\lambda')
        grid on
    end
end
   
if ~PLOT_GRADIENTS
    % Plot the solution
    N = 100;
    DT = T/((N)*M_plot);
    
    X_traj = zeros(4, N+1);
    X_traj(:,1) = [x01; 0; w_opt(1); w_opt(2);];
    
    time = [0:N]*DT*M_plot;
    for i = 2:N+1
        X_traj(:,i) = full(Fs(X_traj(:,i-1), DT));
    end
    
    plot(X_traj(1,:), X_traj(2,:))
    hold all
    xlabel('y')
    ylabel('z')
    grid on
end