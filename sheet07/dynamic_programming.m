clc;
clear all;
close all;

Ts = 0.1;
nSteps = 1;
Q = diag([100,0.01]);
R = 0.001;
x1_max = 2*pi;
x1_min = -pi/2;
x2_max = 10;
u_max = 10;

PERFORM_LQR_DESIGN = 0;

if PERFORM_LQR_DESIGN
    LQR_design(Q, R, Ts, nSteps);
end

N = 20;
input.Ts = Ts; 
input.nSteps = nSteps; 
input.sens = 0; % no sensitivities are needed here

N_x1 = 200; % number of discretization points for state x1
N_x2 = 40;  % number of discretization points for state x2

x1_values = linspace(x1_min,x1_max,N_x1);
x2_values = linspace(-x2_max,x2_max,N_x2);

% mesh for plotting
[X1, X2] = meshgrid(x1_values, x2_values);

N_u = 20; % number of discretization points for control u
u_values = linspace(-u_max,u_max,N_u);

% THESE VARIABLES ARE GENERATED AND STORED BY RUNNING LQR_design.m
load lqr.mat A B Q R K P

%% initialize the cost-to-go function with the cost matrix associated with 
% LQR controller

LQR_cost = zeros(N_x1, N_x2);
LQR_u = zeros(N_x1, N_x2);

for i = 1:N_x1
    for j = 1:N_x2
        x = [x1_values(i); x2_values(j)];
        % INSERT YOUR CODE HERE:
        LQR_cost(i,j) = x'*P*x;
        LQR_u(i,j) = -K*x;
    end
end

J_cost = LQR_cost;

runtime = [];
for k = N-1:-1:1
    tic
    u_map = NaN*ones(N_x1, N_x2);
    J_new = inf+J_cost;
    for i1 = 1:N_x1
        for i2 = 1:N_x2
            cost_local = zeros(N_u,1);
            x_k = [x1_values(i1); x2_values(i2)];
            input.x = x_k;
            for j = 1:N_u
                u_k = u_values(j);
                input.u = u_k;
                output = RK4_integrator( @ode, input );
                x_next = output.value;
                
                % project on discretization grid
                % INSERT YOUR CODE HERE:
                i1_next = project(x_next(1),x1_values);
                i2_next = project(x_next(2),x2_values);
                
                if i1_next <= 0 || i1_next > N_x1 || i2_next <= 0 || i2_next > N_x2
                    cost = Inf;
                else
                    % INSERT YOUR CODE HERE:
                    cost = x_k'*Q*x_k+u_k'*R*u_k + J_cost(i1_next,i2_next);
                end
                
                cost_local(j) = cost;
                
                if cost < J_new(i1,i2)
                    % INSERT YOUR CODE HERE:
                    % Equation (8.4)
                    %P_new = Q + A'*P*A - (S'+A'*P*B)*K;
                    J_new(i1,i2) = cost;
                    u_map(i1,i2) = u_k;
                end % end if cost < J_new
            end % end for j
        end % end for i2
    end % end for i1
    disp(sprintf('k = %d\t%f', k, toc));
    
    J_cost = J_new;
    
    figure(1);
    clf;
    subplot(211);
    surf(X1.', X2.', J_cost); hold all;
    plot3(X1.', X2.', LQR_cost,'--r');
    xlabel('x_1'); ylabel('x_2'); zlabel('J_{cost}');
    legend('DP', 'LQR')
    title(['Cost-to-go function for k = ' num2str(k)])
    drawnow
    
    subplot(212);
    surf(X1.', X2.', u_values(u_map); hold on;
    plot3(X1.', X2.', LQR_u,'--r');
    xlabel('x_1'); ylabel('x_2'); zlabel('u_{map}');
    zlim([-u_max; u_max])
    legend('DP', 'LQR')
    title(['Optimal feedback control for k = ' num2str(k)])
    drawnow
    
end

save DP.mat x1_values x2_values u_values J_cost u_map Ts nSteps u_max

