%% clear system and vars
clc; 
clear;
close all;
cvx_clear;

%% user Input
% define the system
input_dim = 2; % size of input vector U_i
state_dim = 4; % size of state vector (position and velocity) X_i

% bounds on input space
umax = 2; % max absolute value for an input

% time steps
time_horizon = 15; % number of steps from initial condition to completion
time_step = 450 / time_horizon; % 6 minutes total split into even sections

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

%%  optimization parameters
% max iterations
kmax = 100;

% collision avoid region radius
r = 10;

% probailities of safety
alpha = .9; % collision avoidance probability
beta = .9; % goal achievement

% convergence perameters
epsilon_dc = 1e-8; % convergence in cost
epsilon_lambda = 1e-8; % convergence of sum of slack variables to zero

% cost of slack variable
tau_max = 100;
gamma = 1.2;
tau = min(tau_max * ones(kmax,1), gamma.^(0:(kmax-1))');

% storage initial cost for convergence check
input_cost = [1e10; zeros(kmax,1)];
lambda_sum = [1e10; zeros(kmax,1)];
total_cost = [1e20; zeros(kmax,1)];

% initial conditions for U_i
U_a_p = ones(time_horizon * input_dim,1);
U_b_p = ones(time_horizon * input_dim,1);
U_c_p = ones(time_horizon * input_dim,1);


%% generate LTI systems
% LTI variables assumed IID
params = CwhSystemParameters('SamplingPeriod', time_step);

sys = getCwhLtiSystem(4, ...
                      Polyhedron('lb', -umax*ones(input_dim,1), ...
                                 'ub',  umax*ones(input_dim,1)), ...
                      RandomVector('Gaussian', ...
                                   zeros(state_dim,1), ...
                                   diag([1e-4, 1e-4, 5e-8, 5e-8])), ...
                      params);

sysnoi = LtvSystem('StateMatrix',sys.state_mat, ...
                   'DisturbanceMatrix', sys.dist_mat, ...
                   'Disturbance',sys.dist);

% polytope representation of \mathcal{U}
[concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys, time_horizon);

% compute the input concatenated transformations
[A, Cu, Cw] = getConcatMats(sys, time_horizon);

% compute mean_X_sans_input, cov_X_sans_input

X_a_sans_input_rv = SReachFwd('concat-stoch', sysnoi, x_0_a, time_horizon);
X_b_sans_input_rv = SReachFwd('concat-stoch', sysnoi, x_0_b, time_horizon);
X_c_sans_input_rv = SReachFwd('concat-stoch', sysnoi, x_0_c, time_horizon);

mean_X_a_sans_input = X_a_sans_input_rv.mean();
mean_X_b_sans_input = X_b_sans_input_rv.mean();
mean_X_c_sans_input = X_c_sans_input_rv.mean();

mean_X_a_sans_input = mean_X_a_sans_input(sysnoi.state_dim+1:end);
mean_X_b_sans_input = mean_X_b_sans_input(sysnoi.state_dim+1:end);
mean_X_c_sans_input = mean_X_c_sans_input(sysnoi.state_dim+1:end);

cov_X_sans_input = X_a_sans_input_rv.cov();
cov_X_sans_input = cov_X_sans_input(sysnoi.state_dim+1:end, sysnoi.state_dim+1:end);

%% compute the number of polytopic halfspaces to worry about
n_lin_state_a = size(target_set_a.A,1);
n_lin_state_b = size(target_set_b.A,1);
n_lin_state_c = size(target_set_c.A,1);


%% compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i} = ||\sqrt\Sigma_X*h_i||
[sqrt_cov_X_sans_input, p] = chol(cov_X_sans_input(end-3:end,end-3:end));

if p > 0
    sqrt_cov_X_sans_input = sqrt(cov_X_sans_input);
end

scaled_sigma_a_vec = norms(target_set_a.A * sqrt_cov_X_sans_input',2,2);
scaled_sigma_b_vec = norms(target_set_b.A * sqrt_cov_X_sans_input',2,2);
scaled_sigma_c_vec = norms(target_set_c.A * sqrt_cov_X_sans_input',2,2);

%% lower bound for collision avoidance constraint
% compute ||A^1/2||_lb s.t. ||A^1/2||_lb ||x|| <= ||Ax||
% ||A^1/2||_lb is the smallest eigenvalue of A^1/2
sigma_norm_lb = zeros(time_horizon, 1);
for i = 1:time_horizon
    index = 4*(i-1) + (1:2); 
    e = eig(2 * cov_X_sans_input(index, index));
    sigma_norm_lb(i) =  min(e);
end

% bound = r + \Phi^{-1}_{Rayl}(\alpha)||(Sigma_a[t] + Sigma_b[t])^1/2||_lb
collision_avoid_lb = r + raylinv(alpha, 1)*sqrt(sigma_norm_lb);
collision_avoid_lb_sq = collision_avoid_lb.^2;

%% obtain the piecewise linear overapproximation of norminvcdf in [0,0.5]
pwa_accuracy = 1e-3; % Set the maximum piecewise-affine overapproximation error to 1e-3
[invcdf_approx_m, invcdf_approx_c, lb_delta_i] = computeNormCdfInvOverApprox(0.5, pwa_accuracy, n_lin_state_a);

%% initial mean vector for given initial condition
mean_X_a = mean_X_a_sans_input + Cu * U_a_p;
mean_X_b = mean_X_b_sans_input + Cu * U_b_p;
mean_X_c = mean_X_c_sans_input + Cu * U_c_p;

%% Set defaults for cvx
cvx_solver Gurobi
cvx_precision best

%% iterative optimization problem
tic;
k = 1;
while k <= kmax 
    fprintf('iteration: %d ', k);

    % update collision avoid probabilities and gradient
    [g_ab, del_g_ab] = update_g(mean_X_a, mean_X_b, Cu, time_horizon);
    [g_ac, del_g_ac] = update_g(mean_X_a, mean_X_c, Cu, time_horizon);
    [g_bc, del_g_bc] = update_g(mean_X_b, mean_X_c, Cu, time_horizon);

    cvx_begin quiet
        variable U_a(sys.input_dim * time_horizon,1);
        variable U_b(sys.input_dim * time_horizon,1);
        variable U_c(sys.input_dim * time_horizon,1);

        variable mean_X_a(sys.state_dim * time_horizon, 1);
        variable mean_X_b(sys.state_dim * time_horizon, 1);
        variable mean_X_c(sys.state_dim * time_horizon, 1);

        variable lambda_i_ab(time_horizon, 1);
        variable lambda_i_ac(time_horizon, 1);
        variable lambda_i_bc(time_horizon, 1);

        variable delta_i_a(n_lin_state_a, 1);
        variable delta_i_b(n_lin_state_b, 1);
        variable delta_i_c(n_lin_state_c, 1);

        variable norminvover_a(n_lin_state_a, 1);
        variable norminvover_b(n_lin_state_b, 1);
        variable norminvover_c(n_lin_state_c, 1);

        minimize (tau(k)*(sum(lambda_i_ab) + sum(lambda_i_ac) + sum(lambda_i_bc)) + U_a'*U_a + U_b'*U_b + U_c'*U_c)
        subject to
            %----------------------------
            % linear equations defining the state
            %----------------------------
            mean_X_a == mean_X_a_sans_input + Cu * U_a;
            mean_X_b == mean_X_b_sans_input + Cu * U_b;
            mean_X_c == mean_X_c_sans_input + Cu * U_c; 

            %----------------------------
            % u \in \mathcal{U} 
            %----------------------------
            concat_input_space_A * U_a <= concat_input_space_b;
            concat_input_space_A * U_b <= concat_input_space_b; 
            concat_input_space_A * U_c <= concat_input_space_b;

            %----------------------------
            % colission avoidance constraint
            %----------------------------

            % difference of convex function representation of 
            % 0 - (- ||x_a - x_b||) >= r + Raylinv(1-\alpha)*||(Sigma_a + Sigma_b)^1/2||_lb - slack
            % slack variables added for feasibility.
            lambda_i_ab >= 0;
            lambda_i_ac >= 0;
            lambda_i_bc >= 0;

            g_ab + del_g_ab * [U_a - U_a_p;U_b - U_b_p] >= ...
                collision_avoid_lb_sq - lambda_i_ab;
            g_ac + del_g_ac * [U_a - U_a_p;U_c - U_c_p] >= ...
                collision_avoid_lb_sq - lambda_i_ac;
            g_bc + del_g_bc * [U_b - U_b_p;U_c - U_c_p] >= ...
                collision_avoid_lb_sq - lambda_i_bc;

            %----------------------------
            % terminal state constraint
            %----------------------------

            % approximation of inverse normal in convex region
            for delta_i_indx = 1:n_lin_state_a
                norminvover_a(delta_i_indx) >= invcdf_approx_m.* delta_i_a(delta_i_indx) + invcdf_approx_c;
            end
            for delta_i_indx = 1:n_lin_state_b
                norminvover_b(delta_i_indx) >= invcdf_approx_m.* delta_i_b(delta_i_indx) + invcdf_approx_c;
            end
            for delta_i_indx = 1:n_lin_state_c
                norminvover_c(delta_i_indx) >= invcdf_approx_m.* delta_i_c(delta_i_indx) + invcdf_approx_c;
            end

            % \mu_v in target shrunk by \beta
            target_set_a.A * mean_X_a(end-3:end) + scaled_sigma_a_vec.* norminvover_a <= target_set_a.b;
            target_set_b.A * mean_X_b(end-3:end) + scaled_sigma_b_vec.* norminvover_b <= target_set_b.b;
            target_set_c.A * mean_X_c(end-3:end) + scaled_sigma_c_vec.* norminvover_c <= target_set_c.b;

            % \delta_i,v not infinity
            delta_i_a >= lb_delta_i;
            delta_i_b >= lb_delta_i;
            delta_i_c >= lb_delta_i;

            % \delta_i,v in convex region
            delta_i_a <= 0.5;
            delta_i_b <= 0.5;
            delta_i_c <= 0.5;

            % Prob(vehicles \in targets) <= 1-\beta
            sum(delta_i_a) <= 1 - beta;
            sum(delta_i_b) <= 1 - beta;
            sum(delta_i_c) <= 1 - beta;
    cvx_end

    % update Costs
    input_cost(k+1) = U_a'*U_a + U_b'*U_b + U_c'*U_c;
    lambda_sum(k+1) = sum(lambda_i_ab) + sum(lambda_i_ac) + sum(lambda_i_bc);
    total_cost(k+1) = cvx_optval;

    % calculate convergence criteria
    conv_check = abs(input_cost(k+1) - input_cost(k) + tau(k)*(lambda_sum(k+1) - lambda_sum(k)));

    % print statistics
    fprintf('\t %e', conv_check);
    fprintf('\t %e \n', lambda_sum(k+1));

    % check for solved status
    if strcmpi(cvx_status, 'Solved')
        % check for convergence
        if (conv_check <= epsilon_dc) && (lambda_sum(k+1) <= epsilon_lambda)                 
           break
        end

        % if not converged update previous answer to current answer
        U_a_p = U_a;
        U_b_p = U_b;
        U_c_p = U_c;

    % if results are NaN break before error
    elseif strcmpi(cvx_status, 'Failed') || strcmpi(cvx_status, 'Infeasible')
        break
    end

    % update itteration number
    k = k + 1;
end

% check time 
time = toc;

% make k not less than or equal to kmax
k = min(k, kmax);

%% print some useful information
fprintf('\n %s \n', cvx_status);
fprintf('Computation time (sec): %f \n', time);
fprintf('Optimal Total Cost: %f \n', total_cost(k+1));
fprintf('Optimal Slack Cost: %f \n', lambda_sum(k+1));
fprintf('Optimal Input Cost: %f \n', input_cost(k+1));

%% verification
if strcmpi(cvx_status, 'Solved')
    problem(1).initial_condition = x_0_a;
    problem(2).initial_condition = x_0_b;
    problem(3).initial_condition = x_0_c;

    problem(1).input = U_a;
    problem(2).input = U_b;
    problem(3).input = U_c;

    problem(1).target_set = target_set_a;
    problem(2).target_set = target_set_b;
    problem(3).target_set = target_set_c;

    verification_dc = verify(10e5, sys, time_horizon, "L2", r, problem);
end


%% run particle control
particle_control

%% make graphs Our method
motion_path_graph( ...
    [x_0_a; mean_X_a], [x_0_b; mean_X_b], [x_0_c; mean_X_c],...
    [repmat(x_0_a,1,N); mean_X_a], [x_0_b; mean_X_b], [x_0_c; mean_X_c],...
    target_set_a, target_set_b, target_set_c, 1);
distance_graph(...
    [x_0_a; mean_X_a], [x_0_b; mean_X_b], [x_0_c; mean_X_c],...
    [x_0_a; mean_X_a], [x_0_b; mean_X_b], [x_0_c; mean_X_c],...
    r, time_horizon, collision_avoid_lb, 1);
cost_graph
cum_cost(...
    U_a, U_b, U_c,...
    U_a_bl, U_b_bl, U_c_bl, ...
    time_horizon);