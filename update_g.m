function [g, eps, del_g] = update_g(mu_1, mu_2, sigma_1, sigma_2, Cu, r)
    mu = mu_1 - mu_2;
    mupr = mu + r;
    mumr = mu - r;
    sigma = sigma_1 + sigma_2;
    sigmad = diag(sigma);
        
    P_x1 = normcdf(mupr(1:4:end) ./ sqrt(sigmad(1:4:end))) - normcdf(mumr(1:4:end) ./ sqrt(sigmad(1:4:end)));
    P_x2 = normcdf(mupr(2:4:end) ./ sqrt(sigmad(2:4:end))) - normcdf(mumr(2:4:end) ./ sqrt(sigmad(2:4:end)));
    p_x1 = normpdf(mupr(1:4:end) ./ sqrt(sigmad(1:4:end))) - normpdf(mumr(1:4:end) ./ sqrt(sigmad(1:4:end)));
    p_x2 = normpdf(mupr(2:4:end) ./ sqrt(sigmad(2:4:end))) - normpdf(mumr(2:4:end) ./ sqrt(sigmad(2:4:end)));
   
    P_in = P_x1.*P_x2;
        
    C_x1 = Cu(1:4:end,:);
    C_x2 = Cu(2:4:end,:);
    
    del_g_1 = (C_x2 ./ sqrt(sigmad(2:4:end)) .* (P_x1 .* -p_x2)) + ...
              (C_x1 ./ sqrt(sigmad(1:4:end)) .* (P_x2 .* -p_x1));
    del_g_2 = (C_x2 ./ sqrt(sigmad(2:4:end)) .* (P_x1 .*  p_x2)) + ...
              (C_x1 ./ sqrt(sigmad(1:4:end)) .* (P_x2 .*  p_x1));

    g = P_in;
    eps = max(P_in, 1e-6*ones(size(P_in)));
    del_g = [del_g_1, del_g_2];
end