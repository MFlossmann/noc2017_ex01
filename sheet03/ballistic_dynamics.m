function [ dx ] = ballistic_dynamics( x ) 

    % Insert your code here 
    w = 2; % side wind velocity
    d = 0.1; % drag
    g = 9.81; % gravitational acceleration
    
    dx = [x(3); ...
          x(4); ...
          -(x(3) - w)* norm(x(3:4) - [w;0])*d; ...
          -x(4)*norm(x(3:4) - [w;0])*d - g];

end

