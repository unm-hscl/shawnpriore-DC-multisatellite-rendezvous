function [log_g, del_log_g] = update_g(mu_1, mu_2, sigma_1, sigma_2, Cu, r, time_horizon)
    % calculate and extract input
    mu = mu_1 - mu_2;
    sigma = sigma_1 + sigma_2;
    
    % memory holders
    P = zeros(time_horizon, 1);
    gradient_P = zeros(time_horizon, time_horizon * 2);
    
    % iterate through time index
    for i = 1:time_horizon
        % get relavent indexes
        index = 4*(i-1) + (1:2);
        
        % calculate L_2 norm of mean
        mu_i = mu(index);
        mu_i_2norm = norm(mu_i); 
        
        % calculate frobenius norm of Sigma^{1/2}
        sigma_i = sigma(index,index);
        sigma_i_chol = chol(sigma_i);
        sigma_i_fnorm = norm(sigma_i_chol, 'fro');
        
        % get indexed row of controlability matrix
        Cu_i_1 = Cu(index(1), :);
        Cu_i_2 = Cu(index(2), :);
        
        % calculate Prob{|| \nu || <= (r - || \mu ||) / || \Sigma^1/2} ||_F }
        P(i) = raylcdf(max(r - mu_i_2norm, 1e-10) / sigma_i_fnorm, 1);

        % calculate gradient of above probability
        gradient_P(i,:) = max(raylpdf((r - mu_i_2norm)/ sigma_i_fnorm, 1), 1e-12) * ...
            (Cu_i_1 * mu_i(1) + Cu_i_2 * mu_i(2)) / ...
            (sigma_i_fnorm * mu_i_2norm);
    end
    
    % log for concave functions
    log_g = log(P);
    del_log_g = [-gradient_P, gradient_P] ./ P;
end

