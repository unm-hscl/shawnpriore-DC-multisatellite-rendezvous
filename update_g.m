function [log_g, del_log_g] = update_g(mu_1, mu_2, sigma_1, sigma_2, Cu, r)
    % calculate and extract input
    mu = mu_1 - mu_2;
    mu_x = mu(1:4:end);
    mu_y = mu(2:4:end);
    
    mu_x_pr = mu_x + r;
    mu_y_pr = mu_y + r;
    mu_x_mr = mu_x - r;
    mu_y_mr = mu_y - r;
    
    sigma = sigma_1 + sigma_2;
    sigmad = diag(sigma);
    sigmad_x = sigmad(1:4:end);
    sigmad_y = sigmad(2:4:end);
    
    min_P = 1e-10;
    min_p = 1e-15;
        
    % Calculate Prob(in collision region on axis)
    P_x = normcdf(mu_x_pr ./ sqrt(sigmad_x)) - normcdf(mu_x_mr ./ sqrt(sigmad_x));
    P_y = normcdf(mu_y_pr ./ sqrt(sigmad_y)) - normcdf(mu_y_mr ./ sqrt(sigmad_y));
    % Calculate derivative of Prob(in collision region on axis)
    p_x = normpdf(mu_x_pr ./ sqrt(sigmad_x)) - normpdf(mu_x_mr ./ sqrt(sigmad_x));
    p_y = normpdf(mu_y_pr ./ sqrt(sigmad_y)) - normpdf(mu_y_mr ./ sqrt(sigmad_y));  
    
    change_x = (p_x < min_p);
    change_y = (p_y < min_p);
	P_x = max(P_x, min_P);
	P_y = max(P_y, min_P);

    if sum(change_x) > 0
        sign_mu_x = -sign(mu_x);
        p_x(change_x) = sign_mu_x(change_x) * min_p;
    end
    if sum(change_y) > 0
        sign_mu_y = -sign(mu_y);
        p_y(change_y) = sign_mu_y(change_y) * min_p;
    end
    
    % Calculate Prob(in collision region)
    P_in = P_x.*P_y;
        
    % Extract line from controllability matrix
    C_x = Cu(1:4:end,:);
    C_y = Cu(2:4:end,:);
    
    % Calculate gradients
    del_g_1 = (C_y ./ sqrt(sigmad(2:4:end)) .* (P_x .* -p_y)) + ...
              (C_x ./ sqrt(sigmad(1:4:end)) .* (P_y .* -p_x)) ./ ...
              P_in;
    del_g_2 = (C_y ./ sqrt(sigmad(2:4:end)) .* (P_x .*  p_y)) + ...
              (C_x ./ sqrt(sigmad(1:4:end)) .* (P_y .*  p_x)) ./ ...
              P_in;
    
    % Return probability and gradient
    log_g = log(P_in);
    del_log_g = [del_g_1, del_g_2];
end

