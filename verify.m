function results = verify(samples, sys, time_horizon, norm_choice, r, problem) 
    % define variables
    x_0_a = problem(1).initial_condition;
    x_0_b = problem(2).initial_condition;
    x_0_c = problem(3).initial_condition;

    U_a = problem(1).input;
    U_b = problem(2).input;
    U_c = problem(3).input;
    
    target_set_a = problem(1).target_set;
    target_set_b = problem(2).target_set;
    target_set_c = problem(3).target_set;

    % get mc samples
    [state_realization_a, ~, ~] = ...
        generateMonteCarloSims(samples, sys, x_0_a, ...
        time_horizon, U_a);

    [state_realization_b, ~, ~] = ...
        generateMonteCarloSims(samples, sys, x_0_b, ...
        time_horizon, U_b);

    [state_realization_c, ~, ~] = ...
        generateMonteCarloSims(samples, sys, x_0_c, ...
        time_horizon, U_c);

    % test goal achievement
    P_goal_a = 0;
    P_goal_b = 0;
    P_goal_c = 0;

    for i = 1:samples
       P_goal_a = P_goal_a + target_set_a.contains( state_realization_a(end-3:end, i) );
       P_goal_b = P_goal_b + target_set_b.contains( state_realization_b(end-3:end, i) );
       P_goal_c = P_goal_c + target_set_c.contains( state_realization_c(end-3:end, i) );
    end

    results.P_goal_a = P_goal_a / samples;
    results.P_goal_b = P_goal_b / samples;
    results.P_goal_c = P_goal_c / samples;
    
    % test collision avoidance
    diff_ab = state_realization_a - state_realization_b;    
    diff_ac = state_realization_a - state_realization_c;
    diff_bc = state_realization_b - state_realization_c;
    
    norms_ab = zeros(time_horizon, samples);
    norms_ac = zeros(time_horizon, samples);
    norms_bc = zeros(time_horizon, samples);
    
    if strcmpi(norm_choice, "L2")
        for t = 1:time_horizon
            for i = 1:samples
                norms_ab(t, i) = norm( diff_ab(4*t + (1:2), i) );
                norms_ac(t, i) = norm( diff_ac(4*t + (1:2), i) );
                norms_bc(t, i) = norm( diff_bc(4*t + (1:2), i) );
            end
        end
    elseif strcmpi(norm_choice, "Linf")
        for t = 1:time_horizon
            for i = 1:samples
                norms_ab(t, i) = norm( diff_ab(4*t + (1:2), i) , Inf );
                norms_ac(t, i) = norm( diff_ac(4*t + (1:2), i) , Inf );
                norms_bc(t, i) = norm( diff_bc(4*t + (1:2), i) , Inf );
            end
        end
    end
    
    results.P_collision_ab = sum( ( sum( (norms_ab >= r) , 1) == time_horizon)) / samples;
    results.P_collision_ac = sum( ( sum( (norms_ac >= r) , 1) == time_horizon)) / samples;
    results.P_collision_bc = sum( ( sum( (norms_bc >= r) , 1) == time_horizon)) / samples;
    
end