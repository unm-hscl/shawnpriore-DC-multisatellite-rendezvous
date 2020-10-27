% clear system and vars
clc;
clear;
close all;
cvx_clear;

% Iteration setup
runs = 25;

costs = zeros(runs, 1);

for run = 1:runs

    %%%%%%%%%%%%%%%%%%%%%%%%
    % problem setup
    %%%%%%%%%%%%%%%%%%%%%%%%

    % define the system
    input_dim = 2; % size of input vector U_i
    state_dim = 4; % size of state vector (position and velocity) X_i

    % bounds on input space
    umax = 2; % max absolute value for an input

    % time steps
    time_horizon = 15; % number of steps from initial condition to completion
    sampling_period = 450 / time_horizon; % 6 minutes total split into even sections

    % LTI variables
    mu = zeros(state_dim,1); % mean distrubance
    sigma = run * diag([1e-4, 1e-4, 5e-8, 5e-8]); % covariance of disturbance

    % initial states
    % format: x, y, x., y.
    x_0_a = [100; 10.5; 0; 0] ; % satellite A
    x_0_b = [106; 00;   0; 0] ; % satellite B
    x_0_c = [ 94; 00;   0; 0] ; % satellite C

    % target sets
    % format: x, y, x., y.
    target_set_a = Polyhedron('lb', [-15;   -10; -0.01; -0.01], ... 
                              'ub', [-10;    -5;  0.01;  0.01]);  
    target_set_b = Polyhedron('lb', [- 2.5;  15; -0.01; -0.01], ...
                              'ub', [  2.5;  20;  0.01;  0.01]);    
    target_set_c = Polyhedron('lb', [ 10;   -10; -0.01; -0.01], ... 
                              'ub', [ 15;    -5;  0.01;  0.01]);   

    % collision avoid region radius
    r = 10;

    % safety thresholds
    alpha = .9; % collision avoidance
    beta = .9; % goal achievement

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimization parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    % input initilizations
    U_a_init = zeros(input_dim * time_horizon, 1);
    U_b_init = zeros(input_dim * time_horizon, 1);
    U_c_init = zeros(input_dim * time_horizon, 1);

    % max iterations
    kmax = 100;

    % convergence perameters
    epsilon_dc = 1e-8; % convergence in cost
    epsilon_lambda = 1e-8; % convergence of sum of slack variables to zero

    % cost of slack variable
    tau_max = 100;
    gamma = 1.2;
    tau = min(tau_max * ones(kmax,1), gamma.^(0:(kmax-1))');

    % Set defaults for cvx
    cvx_solver gurobi
    cvx_precision best

    %% generate LTI systems
    % LTI variables assumed IID
    params = CwhSystemParameters('SamplingPeriod', sampling_period);

    sys = getCwhLtiSystem(4, ...
                          Polyhedron('lb', -umax*ones(input_dim,1), ...
                                     'ub',  umax*ones(input_dim,1)), ...
                          RandomVector('Gaussian', ...
                                       mu, ...
                                       sigma), ...
                          params);

    sysnoi = LtvSystem('StateMatrix',sys.state_mat, ...
                       'DisturbanceMatrix', sys.dist_mat, ...
                       'Disturbance',sys.dist);

    % polytope representation of \mathcal{U}
    [input_space_A, input_space_b] = getConcatInputSpace(sys, time_horizon);

    % compute the input concatenated transformations
    [A, Cu, Cw] = getConcatMats(sys, time_horizon);

    % compute mean trajectories without input
    X_a_no_input_rv = SReachFwd('concat-stoch', sysnoi, x_0_a, time_horizon);
    X_b_no_input_rv = SReachFwd('concat-stoch', sysnoi, x_0_b, time_horizon);
    X_c_no_input_rv = SReachFwd('concat-stoch', sysnoi, x_0_c, time_horizon);

    mean_X_a_no_input = X_a_no_input_rv.mean();
    mean_X_b_no_input = X_b_no_input_rv.mean();
    mean_X_c_no_input = X_c_no_input_rv.mean();

    mean_X_a_no_input = mean_X_a_no_input(sysnoi.state_dim+1:end);
    mean_X_b_no_input = mean_X_b_no_input(sysnoi.state_dim+1:end);
    mean_X_c_no_input = mean_X_c_no_input(sysnoi.state_dim+1:end);

    % compute covariance of x 
    cov_X_no_input = X_a_no_input_rv.cov();
    cov_X_no_input = cov_X_no_input(sysnoi.state_dim+1:end, sysnoi.state_dim+1:end);

    %% compute the number of polytopic halfspaces to worry about
    n_lin_state_a = size(target_set_a.A,1);
    n_lin_state_b = size(target_set_b.A,1);
    n_lin_state_c = size(target_set_c.A,1);


    %% compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i} = ||\sqrt\Sigma_X*h_i||
    [sqrt_cov_X_no_input, p] = chol(cov_X_no_input(end-3:end,end-3:end));

    if p > 0
        sqrt_cov_X_no_input = sqrt(cov_X_no_input);
    end

    scaled_sigma_a_vec = norms(target_set_a.A * sqrt_cov_X_no_input',2,2);
    scaled_sigma_b_vec = norms(target_set_b.A * sqrt_cov_X_no_input',2,2);
    scaled_sigma_c_vec = norms(target_set_c.A * sqrt_cov_X_no_input',2,2);

    %% lower bound for collision avoidance constraint
    % compute ||A^1/2||_lb s.t. ||A^1/2||_lb ||x|| <= ||Ax||
    % ||A^1/2||_lb is the smallest eigenvalue of A^1/2
    sigma_norm_lb = zeros(time_horizon, 1);
    for i = 1:time_horizon
        index = 4*(i-1) + (1:2); 
        e = eig(2 * cov_X_no_input(index, index));
        sigma_norm_lb(i) =  sqrt(min(e));
    end

    % bound = r + \Phi^{-1}_{Rayl}(\alpha)||(Sigma_a[t] + Sigma_b[t])^1/2||_lb
    collision_avoid_lb = r + raylinv(alpha, 1)*sqrt(sigma_norm_lb);
    collision_avoid_lb_sq = collision_avoid_lb.^2;

    %% obtain the piecewise linear overapproximation of norminvcdf in [0,0.5]
    pwa_accuracy = 1e-3; % Set the maximum piecewise-affine overapproximation error to 1e-3
    [invcdf_approx_m, invcdf_approx_c, lb_delta] = computeNormCdfInvOverApprox(0.5, pwa_accuracy, n_lin_state_a);

    %% set up and solve problem

    problem_dc.time_horizon = time_horizon;

    problem_dc.Cu = Cu;

    problem_dc.collision_avoid_lb_sq = collision_avoid_lb_sq;

    problem_dc.tau = tau;
    problem_dc.beta = beta;
    problem_dc.epsilon_dc = epsilon_dc;
    problem_dc.epsilon_lambda = epsilon_lambda;

    problem_dc.input_space_A = input_space_A;
    problem_dc.input_space_b = input_space_b;

    problem_dc.kmax = kmax;

    problem_dc.input_dim = input_dim;
    problem_dc.state_dim = state_dim;

    problem_dc.invcdf_approx_m = invcdf_approx_m;
    problem_dc.invcdf_approx_c = invcdf_approx_c;
    problem_dc.lb_delta = lb_delta;

    problem_dc.mean_X_a_no_input = mean_X_a_no_input;
    problem_dc.mean_X_b_no_input = mean_X_b_no_input;
    problem_dc.mean_X_c_no_input = mean_X_c_no_input;

    problem_dc.n_lin_state_a = n_lin_state_a;
    problem_dc.n_lin_state_b = n_lin_state_b;
    problem_dc.n_lin_state_c = n_lin_state_c;

    problem_dc.scaled_sigma_a_vec = scaled_sigma_a_vec;
    problem_dc.scaled_sigma_b_vec = scaled_sigma_b_vec;
    problem_dc.scaled_sigma_c_vec = scaled_sigma_c_vec;

    problem_dc.target_set_a = target_set_a;
    problem_dc.target_set_b = target_set_b;
    problem_dc.target_set_c = target_set_c;

    problem_dc.U_a_init = U_a_init;
    problem_dc.U_b_init = U_b_init;
    problem_dc.U_c_init = U_c_init;

    results_dc = dc_solver(problem_dc, 1);
    if strcmpi(results_dc.status, 'Solved')
        costs(run) = results_dc.total_cost(end);
    else
        run = run - 1;
        break
    end
end

%% Covariance scaleing

fig = figure();
fig.Units    = 'inches';
fig.Position = [0,1,10,9];

hold on
plot((1:run), costs(1:run), 'b-o', 'LineWidth', 1);
xlabel('Covariance Multiplier');
ylabel('Cost');
axis([1 run 13 15])
drawnow()
hold off

set(gca, 'OuterPosition', [0.01, 0.01, 0.98, 0.93]);

l = legend({'Cost of Controller'},...
    'Orientation','horizontal');

set(l,'Position', [0.37,0.94,0.26,0.05],'Units', 'normalized');

set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     20);
