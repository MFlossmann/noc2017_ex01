function [ P_t, hist ] = ballistic_dynamics_RK4( v_start)

    T = 0.5;
    M = 100;
    DT = T/M;

    % X0: [p_1y, p_1z, v_1y, v1_z, p_2y, p_2z, v_2y, v_2z]
    X0 = [0, 0, v_start(1), v_start(2), 10, 0, v_start(3), v_start(4)];


    % RK4 integrator
    Xf = X0;
    
    %if ~strcmp(class(Xf),'casadi.SX') % if not casadi stuff: record sim!
    %    hist = zeros(M,4);
    %end
    for j=1:M
        %if ~strcmp(class(Xf), 'casadi.SX')
        %    hist(1,:) = Xf.to_double();
        %end        
        k_1 = ode(Xf);
        k_2 = ode(Xf + k_1'*DT/2);
        k_3 = ode(Xf + k_2'*DT/2);
        k_4 = ode(Xf + k_3'*DT);
        
        Xf = Xf + DT * (k_1' + 2*k_2' + 2*k_3' + k_4') / 6;
    end

    P_t = [Xf(1), Xf(2), Xf(5), Xf(6)];

end

