function results = particle_control(problem, N)
    % Blackmore 2011 for goal achievement and collision avoidance
    % Primary Coder: Vignesh Sivaramakrishnan
    % Modified by: Shawn Priore

    time_horizon = problem.time_horizon;
    r = problem.r;
    
    x_0_a = problem.x_0_a;
    x_0_b = problem.x_0_b;
    x_0_c = problem.x_0_c;
    
    target_set_a = problem.target_set_a;
    target_set_b = problem.target_set_b;
    target_set_c = problem.target_set_c;
    
    input_dim = problem.input_dim;
    state_dim = problem.state_dim;
    
    dist_mean = problem.mean;
    dist_cov = problem.cov;
    
    alpha = problem.alpha;
    beta = problem.beta;
    
    A = problem.A;
    Cu = problem.Cu;
    Cw = problem.Cw;
    
    input_space_A = problem.input_space_A;
    input_space_b = problem.input_space_b;
    
    % parameters
    % big M arbitrary constant
    large_constant = 5000;

    % polytoupe defining ||x_i-x_j||_inf = r
    Avoid_A = ...
        [1, 0;...
        -1, 0; ...
         0, 1;...
         0, -1];
    Avoid_b = [r; r; r; r];

    % randomly generate the disturbance vector from the standard normal.
    wvec = RandomVector('Gaussian', kron(ones(time_horizon,1), dist_mean), kron(eye(time_horizon), dist_cov)); 
    Wa = wvec.getRealizations(N);
    Wb = wvec.getRealizations(N);
    Wc = wvec.getRealizations(N);

    %% Run optimization problem for an optimal control policy
    % We run an optimization problem to determine the control policy over the
    % time horizon T.

    tic;
    cvx_begin quiet
        variable U_a(input_dim * time_horizon,1);
        variable U_b(input_dim * time_horizon,1);
        variable U_c(input_dim * time_horizon,1);

        variable mean_X_a(state_dim * time_horizon, 1);
        variable mean_X_b(state_dim * time_horizon, 1);
        variable mean_X_c(state_dim * time_horizon, 1);    

        variable x_a(state_dim * time_horizon,N);
        variable x_b(state_dim * time_horizon,N);
        variable x_c(state_dim * time_horizon,N);

        variable t_a(N) binary;
        variable t_b(N) binary;
        variable t_c(N) binary;

        variable c_ab(time_horizon * size(Avoid_A,1), N) binary;
        variable c_ac(time_horizon * size(Avoid_A,1), N) binary;
        variable c_bc(time_horizon * size(Avoid_A,1), N) binary;

        variable sc_ab(time_horizon , N) binary;
        variable sc_ac(time_horizon , N) binary;
        variable sc_bc(time_horizon , N) binary;

        variable ssc_ab(N) binary;
        variable ssc_ac(N) binary;
        variable ssc_bc(N) binary;

        minimize (U_a'*U_a + U_b'*U_b + U_c'*U_c);

        subject to
            mean_X_a == A * x_0_a + Cu * U_a;
            mean_X_b == A * x_0_b + Cu * U_b;
            mean_X_c == A * x_0_c + Cu * U_c;

            x_a(:,1:N) == Cw * Wa + repmat(mean_X_a,1,N);
            x_b(:,1:N) == Cw * Wb + repmat(mean_X_b,1,N);
            x_c(:,1:N) == Cw * Wc + repmat(mean_X_c,1,N);

            input_space_A * U_a <= input_space_b;
            input_space_A * U_b <= input_space_b; 
            input_space_A * U_c <= input_space_b;

            for i = 1:N

                target_set_a.A * x_a(end-3:end,i) - target_set_a.b <= large_constant*t_a(i);
                target_set_b.A * x_b(end-3:end,i) - target_set_b.b <= large_constant*t_b(i);
                target_set_c.A * x_c(end-3:end,i) - target_set_c.b <= large_constant*t_c(i);

                for t = 1:time_horizon

                    Avoid_A * (x_a(4*(t-1) + (1:2),i) - x_b(4*(t-1) + (1:2),i)) - Avoid_b + large_constant * c_ab(4*(t-1) + (1:4) , i) >= 0; 
                    Avoid_A * (x_a(4*(t-1) + (1:2),i) - x_c(4*(t-1) + (1:2),i)) - Avoid_b + large_constant * c_ac(4*(t-1) + (1:4) , i) >= 0; 
                    Avoid_A * (x_b(4*(t-1) + (1:2),i) - x_c(4*(t-1) + (1:2),i)) - Avoid_b + large_constant * c_bc(4*(t-1) + (1:4) , i) >= 0; 

                    sum(c_ab(4*(t-1) + (1:4) , i)) - (size(Avoid_A,1) - 1) <= sc_ab(t, i);
                    sum(c_ac(4*(t-1) + (1:4) , i)) - (size(Avoid_A,1) - 1) <= sc_ac(t, i);
                    sum(c_bc(4*(t-1) + (1:4) , i)) - (size(Avoid_A,1) - 1) <= sc_bc(t, i);
                end

                sum(sc_ab(:,i)) <= time_horizon * ssc_ab(i);
                sum(sc_ac(:,i)) <= time_horizon * ssc_ac(i);
                sum(sc_bc(:,i)) <= time_horizon * ssc_bc(i);
            end

            1/N * sum(t_a) <= 1 - alpha;
            1/N * sum(t_b) <= 1 - alpha;
            1/N * sum(t_c) <= 1 - alpha;

            1/N * sum(ssc_ab) <= 1-beta;
            1/N * sum(ssc_ac) <= 1-beta;
            1/N * sum(ssc_bc) <= 1-beta;

    cvx_end;
    time = toc;
    
    % generate results
    if strcmpi(cvx_status, 'Solved')
        results.status = 'Solved';
        
        results.completion_time = time;
        
        results.U_a = U_a;
        results.U_b = U_b;
        results.U_c = U_c;
        
        results.x_a = x_a;
        results.x_b = x_b;
        results.x_c = x_c;
	elseif strcmpi(cvx_status, 'Inaccurate/Solved')
        results.status = 'Could not find exact solution';
    elseif strcmpi(cvx_status, 'Failed') || strcmpi(cvx_status, 'Infeasible')
        results.status = 'Infeasible';
    end
end
