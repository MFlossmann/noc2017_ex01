function [ P_t ] = ballistic_sim( v_start, T )
    M = 100;
    DT = T/M;

    % X0: [p_1y, p_1z, v_1y, v1_z, p_2y, p_2z, v_2y, v_2z]
    X0 = [0; 0; v_start(1); v_start(2)];

    % RK4 integrator
    Xf = X0;
    
    %Xf = ode45(@ode(Xf), [0,T], X0);
    
     for j=1:M
         k_1 = ballistic_dynamics(Xf);
         k_2 = ballistic_dynamics(Xf + k_1*DT/2);
         k_3 = ballistic_dynamics(Xf + k_2*DT/2);
         k_4 = ballistic_dynamics(Xf + k_3*DT);
         
         Xf = Xf + DT * (k_1 + 2*k_2 + 2*k_3 + k_4) / 6;
     end

    P_t = [Xf(1), Xf(2)];

end

