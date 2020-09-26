%% Blackmore TRo 2011 Code to stay in a feasible set. 
% Primary Coder: Vignesh Sivaramakrishnan
% Modified by: Shawn Priore

%% parameters
% number of particles
N = 10; 

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
wvec = RandomVector('Gaussian', zeros(state_dim * time_horizon, 1), kron(eye(time_horizon), sys.dist.cov)); 
Wa = wvec.getRealizations(N);
Wb = wvec.getRealizations(N);
Wc = wvec.getRealizations(N);

%% Run optimization problem for an optimal control policy
% We run an optimization problem to determine the control policy over the
% time horizon T.
    
tstart = tic;

cvx_begin 
    variable U_a_bl(sys.input_dim * time_horizon,1);
    variable U_b_bl(sys.input_dim * time_horizon,1);
    variable U_c_bl(sys.input_dim * time_horizon,1);

    variable mean_X_a_bl(sys.state_dim * time_horizon, 1);
    variable mean_X_b_bl(sys.state_dim * time_horizon, 1);
    variable mean_X_c_bl(sys.state_dim * time_horizon, 1);    

    variable x_a_bl(sys.state_dim * time_horizon,N);
    variable x_b_bl(sys.state_dim * time_horizon,N);
    variable x_c_bl(sys.state_dim * time_horizon,N);

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

    minimize (U_a_bl'*U_a_bl + U_b_bl'*U_b_bl + U_c_bl'*U_c_bl);

    subject to
        mean_X_a_bl == A * x_0_a + Cu * U_a_bl;
        mean_X_b_bl == A * x_0_b + Cu * U_b_bl;
        mean_X_c_bl == A * x_0_c + Cu * U_c_bl;

        x_a_bl(1:end,1:N) == Cw * Wa + repmat(mean_X_a_bl,1,N);
        x_b_bl(1:end,1:N) == Cw * Wb + repmat(mean_X_b_bl,1,N);
        x_c_bl(1:end,1:N) == Cw * Wc + repmat(mean_X_c_bl,1,N);

        concat_input_space_A * U_a_bl <= concat_input_space_b;
        concat_input_space_A * U_b_bl <= concat_input_space_b; 
        concat_input_space_A * U_c_bl <= concat_input_space_b;

        for i = 1:N

            target_set_a.A * x_a_bl(end-3:end,i) - target_set_a.b <= large_constant*t_a(i);
            target_set_b.A * x_b_bl(end-3:end,i) - target_set_b.b <= large_constant*t_b(i);
            target_set_c.A * x_c_bl(end-3:end,i) - target_set_c.b <= large_constant*t_c(i);

            for t = 1:time_horizon

                Avoid_A * (x_a_bl(4*(t-1) + (1:2),i) - x_b_bl(4*(t-1) + (1:2),i)) - Avoid_b + large_constant * c_ab(4*(t-1) + (1:4) , i) >= 0; 
                Avoid_A * (x_a_bl(4*(t-1) + (1:2),i) - x_c_bl(4*(t-1) + (1:2),i)) - Avoid_b + large_constant * c_ac(4*(t-1) + (1:4) , i) >= 0; 
                Avoid_A * (x_b_bl(4*(t-1) + (1:2),i) - x_c_bl(4*(t-1) + (1:2),i)) - Avoid_b + large_constant * c_bc(4*(t-1) + (1:4) , i) >= 0; 

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
blackmore_time_to_solve = toc(tstart);

%% verification
if strcmpi(cvx_status,'Solved')
    problem(1).input = U_a_bl;
    problem(2).input = U_b_bl;
    problem(3).input = U_c_bl;

    verification_bl = verify(10e5, sys, time_horizon, "linf", r, problem);
end
