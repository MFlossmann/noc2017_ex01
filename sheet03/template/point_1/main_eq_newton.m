%% Implementation of a Gauss-Neston SQP solver usign CasADi

close all
clear variables
clc

HESSIAN_APPROXIMATION = 'exact';
% HESSIAN_APPROXIMATION = 'identity'; rho = 100;

% import CasADi
import casadi.*

nv = 2;
x = MX.sym('x',nv);

% Insert your code here (define objective function and its gradient and Hessian)
% ...
% ...
% ...

% Insert your code here (define constraints and its Jacobian and Hessian)
% ...
% ...
% ...

% SQP solver
max_it = 100;
iter = zeros(nv+1,max_it);
iter(:,1) = [1 -1 1]'; % Initial guess

for i=2:max_it    
    
    x_k = iter(1:2,i-1);
    lambda_k = iter(3,i-1);
    switch HESSIAN_APPROXIMATION
        case 'exact'
        % Insert your code here (exact Hessian)
        % ...
        % ...
        case 'identity'
        % Insert your code here (scaled identity approximation)
        % ...
        % ... 
    end

    % Build and solve the KKT system
    % ...
    % ...
        
end

iter = iter(:,1:i-1);
[X,Y] = meshgrid(-1.5:.05:1.5, -1.5:.05:1.5);
Z = log(1 + 1/2*(X -1).^2 + 1/2*(10*(Y -X.^2)).^2 + 1/2*Y.^2);
figure()
subplot(1,2,1)
surf(X,Y,Z)
xlabel('x_1')
ylabel('x_2')
hold all
plot3(iter(1,:),iter(2,:),zeros(length(iter(1,:)),1),'black')
plot3(iter(1,:),iter(2,:),zeros(length(iter(1,:)),1),'blacko')
y_g = linspace(-0.08,1.1,20);
x_g = -(1 - y_g).^2;
plot3(x_g,y_g,zeros(length(x_g),1),'r');

x_h = linspace(-1.5,1.5,20);
y_h = x_h.^2;
plot3(x_h,y_h,zeros(length(x_h),1),'r--');

xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
subplot(1,2,2)
contour(X,Y,Z)
hold all
plot3(iter(1,:),iter(2,:),zeros(length(iter(1,:)),1),'black')
plot3(iter(1,:),iter(2,:),zeros(length(iter(1,:)),1),'blacko')
xlabel('x_1')
ylabel('x_2')
plot3(x_g,y_g,zeros(length(x_g),1),'r');
plot3(x_h,y_h,zeros(length(x_h),1),'r--');
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
figure()
plot(iter(1:2,:)')
grid on
xlabel('iterations')
ylabel('primal solution')
grid on