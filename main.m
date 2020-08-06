%% Clear system and vars
clc; 
clearvars;
clear all;

%% User Input
% Define the system
input_dim = 2; % Size of Input Vector U_i
state_dim = 4; % Size of State Vector (position and velocity) X_i

% Bounds on Inpux and State Space Bounds
umax = 1; % Max value for an input

% Time Steps
time_horizon = 20;
time_step = 600 / time_horizon;

% Initial States
% Format: x, y, x., y.
x_0_a = [- 95; -15; 0; 0]; 
x_0_b = [-150;   5; 0; 0]; 
x_0_c = [- 55;  50; 0; 0]; 
   
% Target Sets
% Format: x, y, x., y.

% MP 1: A, B, C 
% MP 2: A, C, B 

target_set_a = Polyhedron('lb', [- 2.5;  20; -0.01; -0.01], ...
                          'ub', [  2.5;  25;  0.01;  0.01]);    
target_set_c = Polyhedron('lb', [-15; -15; -0.01; -0.01], ... 
                          'ub', [-10; -10;  0.01;  0.01]);   
target_set_b = Polyhedron('lb', [ 10; -15; -0.01; -0.01], ... 
                          'ub', [ 15; -10;  0.01;  0.01]);  
 

%% Generate LTI Systems
% LTI Variables assumed IID
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

% Polytope representation of \mathcal{U}
[concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys, time_horizon);

% Compute the input concatenated transformations
[A, Cu, ~] = getConcatMats(sys, time_horizon);
                      
% Compute mean_X_sans_input, cov_X_sans_input

X_a_sans_input_rv = SReachFwd('concat-stoch', sysnoi, x_0_a, time_horizon);
X_b_sans_input_rv = SReachFwd('concat-stoch', sysnoi, x_0_b, time_horizon);
X_c_sans_input_rv = SReachFwd('concat-stoch', sysnoi, x_0_c, time_horizon);

mean_X_a_sans_input = X_a_sans_input_rv.mean();
mean_X_b_sans_input = X_b_sans_input_rv.mean();
mean_X_c_sans_input = X_c_sans_input_rv.mean();

mean_X_a_sans_input = mean_X_a_sans_input(sysnoi.state_dim+1:end);
mean_X_b_sans_input = mean_X_b_sans_input(sysnoi.state_dim+1:end);
mean_X_c_sans_input = mean_X_c_sans_input(sysnoi.state_dim+1:end);

cov_X_a_sans_input = X_a_sans_input_rv.cov();
cov_X_b_sans_input = X_b_sans_input_rv.cov();
cov_X_c_sans_input = X_c_sans_input_rv.cov();

cov_X_a_sans_input = cov_X_a_sans_input(sysnoi.state_dim+1:end, sysnoi.state_dim+1:end);
cov_X_b_sans_input = cov_X_b_sans_input(sysnoi.state_dim+1:end, sysnoi.state_dim+1:end);
cov_X_c_sans_input = cov_X_c_sans_input(sysnoi.state_dim+1:end, sysnoi.state_dim+1:end);

%% Compute M --- the number of polytopic halfspaces to worry about
n_lin_state_a = size(target_set_a.A,1);
n_lin_state_b = size(target_set_b.A,1);
n_lin_state_c = size(target_set_c.A,1);


%% Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i} = ||\sqrt\Sigma_X*h_i||
[sqrt_cov_X_a_sans_input, p_a] = chol(cov_X_a_sans_input(end-3:end,end-3:end));
[sqrt_cov_X_b_sans_input, p_b] = chol(cov_X_b_sans_input(end-3:end,end-3:end));
[sqrt_cov_X_c_sans_input, p_c] = chol(cov_X_c_sans_input(end-3:end,end-3:end));

if p_a > 0
    sqrt_cov_X_a_sans_input = sqrt(cov_X_a_sans_input);
end
if p_b > 0
    sqrt_cov_X_b_sans_input = sqrt(cov_X_b_sans_input);
end
if p_c > 0
    sqrt_cov_X_c_sans_input = sqrt(cov_X_c_sans_input);
end

scaled_sigma_a_vec = norms(target_set_a.A * sqrt_cov_X_a_sans_input',2,2);
scaled_sigma_b_vec = norms(target_set_b.A * sqrt_cov_X_b_sans_input',2,2);
scaled_sigma_c_vec = norms(target_set_c.A * sqrt_cov_X_c_sans_input',2,2);


%% Obtain the piecewise linear overapproximation of norminvcdf in [0,0.5]
pwa_accuracy = 1e-3; % Set the maximum piecewise-affine overapproximation error to 1e-3
[invcdf_approx_m, invcdf_approx_c, lb_delta_i] = computeNormCdfInvOverApprox(0.5, pwa_accuracy, n_lin_state_a);

%%  Optimization Parameters
kmax = 200;
k = 1;

epsilon_collision = .001;
epsilon_tau = .0001;
epsilon_target = .1;
epsilon_dc = .0001;

tau = [.1; zeros(kmax -1, 1)];
tau_max = 20;
gamma = 1.05;

input_cost = [1e20; zeros(kmax,1)];
lambda_sum = [1e20; zeros(kmax,1)];
total_cost = [1e40; zeros(kmax,1)];
r = 5;

%% Initial Condition
U_a_p = -ones(size(Cu,2),1);
U_b_p = -ones(size(Cu,2),1);
U_c_p = -ones(size(Cu,2),1);
mean_X_a = mean_X_a_sans_input + Cu * U_a_p;
mean_X_b = mean_X_b_sans_input + Cu * U_b_p;
mean_X_c = mean_X_c_sans_input + Cu * U_c_p;

%% Collection of solutions
U_a_storage = zeros(size(Cu,2),kmax);
U_b_storage = zeros(size(Cu,2),kmax);
U_c_storage = zeros(size(Cu,2),kmax);

%% Optimization Problem
cvx_solver mosek
cvx_precision default
tic;
while k <= kmax 
    fprintf('itteration: %d \n', k);
    
    [g_ab, eps_ab, del_g_ab] = update_g(mean_X_a, mean_X_b, cov_X_a_sans_input, cov_X_b_sans_input, Cu, r);
    [g_ac, eps_ac, del_g_ac] = update_g(mean_X_a, mean_X_c, cov_X_a_sans_input, cov_X_c_sans_input, Cu, r);
    [g_bc, eps_bc, del_g_bc] = update_g(mean_X_b, mean_X_c, cov_X_b_sans_input, cov_X_c_sans_input, Cu, r);
    
    
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

        minimize (U_a'*U_a + U_b'*U_b + U_c'*U_c + tau(k)*(sum(lambda_i_ab) + sum(lambda_i_ac) + sum(lambda_i_bc)) )
        subject to
            % Linear equations defining the state
            mean_X_a == mean_X_a_sans_input + Cu * U_a;
            mean_X_b == mean_X_b_sans_input + Cu * U_b;
            mean_X_c == mean_X_c_sans_input + Cu * U_c; 
            
            % u \in \mathcal{U} 
            concat_input_space_A * U_a <= concat_input_space_b;
            concat_input_space_A * U_b <= concat_input_space_b; 
            concat_input_space_A * U_c <= concat_input_space_b;
            
            % Colission avoidance constraint
            lambda_i_ab >= 0;
            lambda_i_ac >= 0;
            lambda_i_bc >= 0;
            
            0 - (-del_g_ab * [U_a - U_a_p;U_b - U_b_p]) <= log(epsilon_collision ./ eps_ab) .* g_ab + lambda_i_ab;
            0 - (-del_g_ac * [U_a - U_a_p;U_c - U_c_p]) <= log(epsilon_collision ./ eps_ac) .* g_ac + lambda_i_ac;
            0 - (-del_g_bc * [U_b - U_b_p;U_c - U_c_p]) <= log(epsilon_collision ./ eps_bc) .* g_bc + lambda_i_bc;
            
            % Terminal state constraint
            for delta_i_indx = 1:n_lin_state_a
                norminvover_a(delta_i_indx) >= invcdf_approx_m.* delta_i_a(delta_i_indx) + invcdf_approx_c;
            end
            for delta_i_indx = 1:n_lin_state_b
                norminvover_b(delta_i_indx) >= invcdf_approx_m.* delta_i_b(delta_i_indx) + invcdf_approx_c;
            end
            for delta_i_indx = 1:n_lin_state_c
                norminvover_c(delta_i_indx) >= invcdf_approx_m.* delta_i_c(delta_i_indx) + invcdf_approx_c;
            end
            
            target_set_a.A * mean_X_a(end-3:end) + scaled_sigma_a_vec.* norminvover_a <= target_set_a.b;
            target_set_b.A * mean_X_b(end-3:end) + scaled_sigma_b_vec.* norminvover_b <= target_set_b.b;
            target_set_c.A * mean_X_c(end-3:end) + scaled_sigma_c_vec.* norminvover_c <= target_set_c.b;

            delta_i_a >= lb_delta_i;
            delta_i_b >= lb_delta_i;
            delta_i_c >= lb_delta_i;
            
            delta_i_a <= 0.5;
            delta_i_b <= 0.5;
            delta_i_c <= 0.5;
            
            sum(delta_i_a) + sum(delta_i_b) + sum(delta_i_c) <= epsilon_target;
    cvx_end
        
    input_cost(k+1) = U_a'*U_a + U_b'*U_b + U_c'*U_c;
    lambda_sum(k+1) = sum(lambda_i_ab) + sum(lambda_i_ac) + sum(lambda_i_bc);
    total_cost(k+1) = cvx_optval;
    
    conv_check = abs(input_cost(k+1) - input_cost(k) + tau(k)*(lambda_sum(k+1) - lambda_sum(k)));
    
    if strcmpi(cvx_status, 'Solved')
        if (conv_check <= epsilon_dc) && (lambda_sum(k+1) <= epsilon_tau)                 
           break
        end
        U_a_p = U_a;
        U_b_p = U_b;
        U_c_p = U_c;
        U_a_storage(:, k) = U_a;
        U_b_storage(:, k) = U_b;
        U_c_storage(:, k) = U_c;
    elseif strcmpi(cvx_status, 'Failed') || strcmpi(cvx_status, 'Infeasible')
        break
    end
    tau(k+1) = min(gamma*tau(k), tau_max);    
    k = k + 1;
end


time = toc;
k = min(k, kmax);

fprintf('Complete \n');
fprintf('Computation time (min): %f \n', time / 60);
fprintf('Locally Optimal Cost: %f \n', total_cost(k+1));

%% Choose itteration to plot (if needed)
% sn = 200;
% mean_X_a = mean_X_a_sans_input + Cu * U_a_storage(:,sn);
% mean_X_b = mean_X_b_sans_input + Cu * U_b_storage(:,sn);
% mean_X_c = mean_X_c_sans_input + Cu * U_c_storage(:,sn); 

%% Make Grpahs
motion_path_graph
cost_graph
distance_graph