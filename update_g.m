function [g, del_g] = update_g(mu_1, mu_2, Cu, time_horizon)
    % calculate and extract input
    mu = mu_1 - mu_2;
    
    % memory holders
    g = zeros(time_horizon, 1);
    gradient_g = zeros(time_horizon, time_horizon * 2);
    
    % iterate through time index
    for i = 1:time_horizon
        % get relavent indexes
        index = 4*(i-1) + (1:2);
        
        % calculate L_2 norm of mean
        mu_i = mu(index);
        g(i) = norm(mu_i); 
        
        % get indexed rows of controlability matrix
        Cu_i = Cu(index, :);
        
        % calculate gradient of norm
        gradient_g(i,:) = mu_i' * Cu_i ./ g(i);
    end
    
    % compile gradient w.r.t u_1 and u_2
    del_g = [gradient_g, -gradient_g];
end

