clc;
clear all;
close all;

clip_LQR = 0;

input.sens = 0; % no sensitivities are needed here
x0 = [pi; 0];

% THESE VARIABLES ARE GENERATED AND STORED BY RUNNING LQR_design.m and 
% dynamic programming.m

load lqr.mat A B Q R K P
load DP.mat x1_values x2_values u_values J_cost u_map Ts nSteps u_max

over_sampling = 0.5;
input.Ts = over_sampling*Ts; input.nSteps = nSteps;
% u_map = imgaussfilt(u_map,1);

% Closed-loop simulation
time = 0;
Tf = 5;
state_LQR = x0;
state_DP = x0;
us_DP = [];
us_LQR = [];
cost_DP = [0];
cost_LQR = [0];
iter = 0;
figure(1);
u_DP_old = 0;

while time(end) < Tf
    
    % optimal feedback law from LQR
    % INSERT YOUR CODE HERE:
    u_LQR = -K*state_LQR(:,end);
    if (clip_LQR ==1)
        min(10,max(-10,u_LQR));
    end
    us_LQR = [us_LQR u_LQR];
    
    cost_LQR = [cost_LQR cost_LQR(end)+u_LQR.'*R*u_LQR+state_LQR(:,end).'*Q*state_LQR(:,end)];
    
    % apply LQR control to nonlinear system
    input.x = state_LQR(:,end);
    input.u = u_LQR;
    output = RK4_integrator( @ode, input );
    state_LQR(:,end+1) = output.value;
    
    % optimal feedback law from DP   
    % INSERT YOU CODE HERE:
    i1 = project(state_DP(1,end), x1_values); 
    i2 = project(state_DP(2,end), x2_values);
    u_DP = u_values(u_map(i1,i2));
        
    us_DP = [us_DP u_DP];
    cost_DP = [cost_DP cost_DP(end)+u_DP.'*R*u_DP+state_DP(:,end).'*Q*state_DP(:,end)];
    
    % apply DP control to nonlinear system
    input.x = state_DP(:,end);
    input.u = u_DP;
    output = RK4_integrator( @ode, input );
    state_DP(:,end+1) = output.value;
    
    % next time step and visualize result
    iter = iter+1;
    time(end+1) = iter*over_sampling*Ts;
    
    % plot results
    clf
    subplot(221);
    plot(time,state_LQR,'r'); hold on;
    plot(time,state_DP,'b');
    xlabel('time(s)');
    ylabel('State x');
    legend('LQR', 'DP')
    
    subplot(223);
    stairs(time(1:end-1),us_LQR,'r'); hold on;
    stairs(time(1:end-1),us_DP,'b');
    xlabel('time(s)');
    ylabel('Control u');
    legend('LQR', 'DP')
    
    subplot(2,2,[2 4])
    plot(time,cost_LQR,'--ro'); hold on;
    plot(time,cost_DP,'--bo');
    xlabel('time(s)');
    ylabel('Closed-loop cost');
    legend('LQR', 'DP')
    
    pause(Ts/2);
end

