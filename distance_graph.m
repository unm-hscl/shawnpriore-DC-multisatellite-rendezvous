%% Graph of distance at each time point
function distance_graph(all_a, all_b, all_c, r, time_horizon, collision_avoid_lb)
    % Create empty vector
    dist_ab = zeros(time_horizon,1);
    dist_ac = zeros(time_horizon,1);
    dist_bc = zeros(time_horizon,1);

    % Find distances between vehicles using L_2ty norm
    for i = 0:time_horizon
       dist_ab(i+1) = norm(all_a((1:2)+(i*4))- all_b((1:2)+(i*4)));
       dist_ac(i+1) = norm(all_a((1:2)+(i*4))- all_c((1:2)+(i*4)));
       dist_bc(i+1) = norm(all_b((1:2)+(i*4))- all_c((1:2)+(i*4)));
    end

    % Plot graph
    figure();
    hold on
    plot((0:time_horizon), dist_ab,'r-o');
    plot((0:time_horizon), dist_ac,'g-o');
    plot((0:time_horizon), dist_bc,'b-o');
    plot((0:time_horizon), r*ones(time_horizon+1,1),'k-');
    plot((0:time_horizon), [r; collision_avoid_lb],'k--');
    drawnow()
    hold off

    % Fancy labels
    xlabel('Time Step, $k$', 'interpreter', 'latex');
    ylabel('Distance (in meters)', 'interpreter', 'latex');
    legend({'$||A-B||_2$','$||A-C||_2$','$||B-C||_2$',strcat('R=', num2str(r)), '$R+\Phi^{-1}_{Rayleigh}(\alpha)\|(2\Sigma(k))^{\frac{1}{2}}\|_\ell$'},...
        'interpreter', 'latex',...
        'Location','northwest');

    % Center graph
    axis([-1 time_horizon+1 0 1.1*max([dist_ab;dist_ac;dist_bc])])
end