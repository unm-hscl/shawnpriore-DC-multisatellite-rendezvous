function results = dc_solver(problem)
    % assign variables
    time_horizon = problem.time_horizon;
    
    Cu = problem.Cu;
    
    collision_avoid_lb_sq = problem.collision_avoid_lb_sq;

    tau = problem.tau;
    beta = problem.beta;

    
    epsilon_dc = problem.epsilon_dc;
    epsilon_lambda = problem.epsilon_lambda;
    
    input_space_A = problem.input_space_A;
    input_space_b = problem.input_space_b;
    
    kmax = problem.kmax;
    
    input_dim = problem.input_dim;
    state_dim = problem.state_dim;
        
    invcdf_approx_m = problem.invcdf_approx_m;
    invcdf_approx_c = problem.invcdf_approx_c;
    lb_delta = problem.lb_delta;

    mean_X_a_no_input = problem.mean_X_a_no_input;
    mean_X_b_no_input = problem.mean_X_b_no_input;
    mean_X_c_no_input = problem.mean_X_c_no_input;
    
    n_lin_state_a = problem.n_lin_state_a;
    n_lin_state_b = problem.n_lin_state_b;
    n_lin_state_c = problem.n_lin_state_c;

    scaled_sigma_a_vec = problem.scaled_sigma_a_vec;
    scaled_sigma_b_vec = problem.scaled_sigma_b_vec;
    scaled_sigma_c_vec = problem.scaled_sigma_c_vec;
    
    target_set_a = problem.target_set_a;
    target_set_b = problem.target_set_b;
    target_set_c = problem.target_set_c;
    
    U_a_p = problem.U_a_init;
    U_b_p = problem.U_b_init;
    U_c_p = problem.U_c_init;
    
    % storage initial cost for convergence check
    input_cost = [1e10; zeros(kmax,1)];
    lambda_sum = [1e10; zeros(kmax,1)];
    total_cost = [1e20; zeros(kmax,1)];

    % iterative optimization problem
    tic;
    k = 1;
    while k <= kmax 
        fprintf('iteration: %d ', k);

        % update collision avoid probabilities and gradient
        if k ~= 1
            [g_ab, del_g_ab] = update_g(mean_X_a, mean_X_b, Cu, time_horizon);
            [g_ac, del_g_ac] = update_g(mean_X_a, mean_X_c, Cu, time_horizon);
            [g_bc, del_g_bc] = update_g(mean_X_b, mean_X_c, Cu, time_horizon);
        else
            [g_ab, del_g_ab] = update_g(mean_X_a_no_input + Cu * U_a_p, mean_X_b_no_input + Cu * U_b_p, Cu, time_horizon);
            [g_ac, del_g_ac] = update_g(mean_X_a_no_input + Cu * U_a_p, mean_X_c_no_input + Cu * U_c_p, Cu, time_horizon);
            [g_bc, del_g_bc] = update_g(mean_X_b_no_input + Cu * U_b_p, mean_X_c_no_input + Cu * U_c_p, Cu, time_horizon);
        end

        cvx_begin quiet
            variable U_a(input_dim * time_horizon,1);
            variable U_b(input_dim * time_horizon,1);
            variable U_c(input_dim * time_horizon,1);

            variable mean_X_a(state_dim * time_horizon, 1);
            variable mean_X_b(state_dim * time_horizon, 1);
            variable mean_X_c(state_dim * time_horizon, 1);

            variable lambda_ab(time_horizon, 1);
            variable lambda_ac(time_horizon, 1);
            variable lambda_bc(time_horizon, 1);

            variable delta_a(n_lin_state_a, 1);
            variable delta_b(n_lin_state_b, 1);
            variable delta_c(n_lin_state_c, 1);

            variable norminvover_a(n_lin_state_a, 1);
            variable norminvover_b(n_lin_state_b, 1);
            variable norminvover_c(n_lin_state_c, 1);

            minimize (tau(k)*(sum(lambda_ab) + sum(lambda_ac) + sum(lambda_bc)) + U_a'*U_a + U_b'*U_b + U_c'*U_c)
            subject to
                %----------------------------
                % linear equations defining the state
                %----------------------------
                mean_X_a == mean_X_a_no_input + Cu * U_a;
                mean_X_b == mean_X_b_no_input + Cu * U_b;
                mean_X_c == mean_X_c_no_input + Cu * U_c; 

                %----------------------------
                % u \in \mathcal{U} 
                %----------------------------
                input_space_A * U_a <= input_space_b;
                input_space_A * U_b <= input_space_b; 
                input_space_A * U_c <= input_space_b;

                %----------------------------
                % colission avoidance constraint
                %----------------------------

                % difference of convex function representation of 
                % ||x_a - x_b||^2 >= (r + Raylinv(\alpha)*min_eig(Sigma_w(k))^1/2)^2 - slack
                % slack variables added for feasibility.
                lambda_ab >= 0;
                lambda_ac >= 0;
                lambda_bc >= 0;

                0 - -(g_ab + del_g_ab * [U_a - U_a_p;U_b - U_b_p]) >= ...
                    collision_avoid_lb_sq - lambda_ab;
                0 - -(g_ac + del_g_ac * [U_a - U_a_p;U_c - U_c_p]) >= ...
                    collision_avoid_lb_sq - lambda_ac;
                0 - -(g_bc + del_g_bc * [U_b - U_b_p;U_c - U_c_p]) >= ...
                    collision_avoid_lb_sq - lambda_bc;

                %----------------------------
                % terminal state constraint
                %----------------------------

                % approximation of inverse normal in convex region
                for delta_indx = 1:n_lin_state_a
                    norminvover_a(delta_indx) >= invcdf_approx_m.* delta_a(delta_indx) + invcdf_approx_c;
                end
                for delta_indx = 1:n_lin_state_b
                    norminvover_b(delta_indx) >= invcdf_approx_m.* delta_b(delta_indx) + invcdf_approx_c;
                end
                for delta_indx = 1:n_lin_state_c
                    norminvover_c(delta_indx) >= invcdf_approx_m.* delta_c(delta_indx) + invcdf_approx_c;
                end

                % \mu_v in target shrunk by \beta
                target_set_a.A * mean_X_a(end-3:end) + scaled_sigma_a_vec.* norminvover_a <= target_set_a.b;
                target_set_b.A * mean_X_b(end-3:end) + scaled_sigma_b_vec.* norminvover_b <= target_set_b.b;
                target_set_c.A * mean_X_c(end-3:end) + scaled_sigma_c_vec.* norminvover_c <= target_set_c.b;

                % \delta_i,v not infinity
                delta_a >= lb_delta;
                delta_b >= lb_delta;
                delta_c >= lb_delta;

                % \delta_i,v in convex region
                delta_a <= 0.5;
                delta_b <= 0.5;
                delta_c <= 0.5;

                % Prob(vehicles \in targets) <= 1-\beta
                sum(delta_a) <= 1 - beta;
                sum(delta_b) <= 1 - beta;
                sum(delta_c) <= 1 - beta;
        cvx_end

        % update Costs
        input_cost(k+1) = U_a'*U_a + U_b'*U_b + U_c'*U_c;
        lambda_sum(k+1) = sum(lambda_ab) + sum(lambda_ac) + sum(lambda_bc);
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

    % print some useful information
    fprintf('%s \n', cvx_status);
    fprintf('Computation time (sec): %f \n', time);
    fprintf('Total Cost: %f \n', total_cost(k+1));
    fprintf('Slack Cost: %f \n', lambda_sum(k+1));
    fprintf('Input Cost: %f \n', input_cost(k+1));
    
    % generate results
    if strcmpi(cvx_status, 'Solved')
        results.status = 'Solved';
        
        results.completion_time = time;
        
        results.U_a = U_a;
        results.U_b = U_b;
        results.U_c = U_c;
        
        results.mean_X_a = mean_X_a;
        results.mean_X_b = mean_X_b;
        results.mean_X_c = mean_X_c;
        
        results.total_cost = total_cost(2:k+1);
        results.lambda_sum = lambda_sum(2:k+1);
        results.input_cost = input_cost(2:k+1);
        
        results.max_iteration = k;
	elseif strcmpi(cvx_status, 'Inaccurate/Solved') || (k == kmax)
        results.status = 'Could not find solution within max iterations';
    elseif strcmpi(cvx_status, 'Failed') || strcmpi(cvx_status, 'Infeasible')
        results.status = 'Infeasible';
    end
end