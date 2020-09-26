%% Graph of distance at each time point
function distance_graph(all_a, all_b, all_c, all_a_pc, all_b_pc, all_c_pc, r, time_horizon, collision_avoid_lb, particles)
    % Create empty vector
    dist_ab = zeros(time_horizon,1);
    dist_ac = zeros(time_horizon,1);
    dist_bc = zeros(time_horizon,1);
    dist_ab_pc = zeros(time_horizon,particles);
    dist_ac_pc = zeros(time_horizon,particles);
    dist_bc_pc = zeros(time_horizon,particles);

    % Find distances between vehicles using L_2ty norm
    for i = 0:time_horizon
       dist_ab(i+1) = norm(all_a((1:2)+(i*4))- all_b((1:2)+(i*4)));
       dist_ac(i+1) = norm(all_a((1:2)+(i*4))- all_c((1:2)+(i*4)));
       dist_bc(i+1) = norm(all_b((1:2)+(i*4))- all_c((1:2)+(i*4)));
       for j = 1:particles
           dist_ab_pc(i+1,j) = norm(all_a_pc((1:2)+(i*4),j)- all_b_pc((1:2)+(i*4),j), 'Inf');
           dist_ac_pc(i+1,j) = norm(all_a_pc((1:2)+(i*4),j)- all_c_pc((1:2)+(i*4),j), 'Inf');
           dist_bc_pc(i+1,j) = norm(all_b_pc((1:2)+(i*4),j)- all_c_pc((1:2)+(i*4),j), 'Inf');
       end
    end
    
    
    pc_mean = zeros(time_horizon+1, 3);
    pc_mean(:,1) = mean(dist_ab_pc, 2);
    pc_mean(:,2) = mean(dist_ac_pc, 2);
    pc_mean(:,3) = mean(dist_bc_pc, 2);
    
    pc_min = zeros(time_horizon+1, 3);
    pc_min(:,1) = min(dist_ab_pc, [], 2);
    pc_min(:,2) = min(dist_ac_pc, [], 2);
    pc_min(:,3) = min(dist_bc_pc, [], 2);
    pc_min = pc_mean - pc_min;
    
    pc_max = zeros(time_horizon+1, 3);
    pc_max(:,1) = max(dist_ab_pc, [], 2);
    pc_max(:,2) = max(dist_ac_pc, [], 2);
    pc_max(:,3) = max(dist_bc_pc, [], 2);
    pc_max = pc_max - pc_mean;
    
    axis_max = max( max(pc_max, [], 'all'), max([dist_ab;dist_ac;dist_bc]));

    % Plot graph
    fig = figure();
    fig.Units    = 'inches';
    fig.Position = [0.75,0,10,11.25];

    
    subplot(2,1,1);
    hold on
    p1 = plot((0:time_horizon), dist_ab, 'Color', [0 0.4470 0.7410], 'Marker', 's');
    p2 = plot((0:time_horizon), dist_ac, 'Color', [0.8500 0.3250 0.0980], 'Marker', 'o');
    p3 = plot((0:time_horizon), dist_bc, 'Color', [0.4940 0.1840 0.5560], 'Marker', 'd');
    p4 = plot((0:time_horizon), r*ones(time_horizon+1,1),'k-', 'LineWidth', 2);
    p5 = plot((0:time_horizon), [r; collision_avoid_lb],'k--');
    drawnow()
    hold off
    xlabel('Time Step, k');
    ylabel('Distance (in meters)');
    axis([-1 time_horizon+1 5 1.1*axis_max])    
    set(gca, 'OuterPosition', [0.01,0.47,0.98,0.45]);
        
    % Plot particle control
    subplot(2,1,2);
    hold on
    errorbar((0:time_horizon), pc_mean(:,1), pc_min(:,1), pc_max(:,1), 'LineStyle', '--', 'Color', [0 0.4470 0.7410], 'Marker', 's');
    errorbar((0:time_horizon), pc_mean(:,2), pc_min(:,2), pc_max(:,2), 'LineStyle', '--', 'Color', [0.8500 0.3250 0.0980], 'Marker', 'o');
    errorbar((0:time_horizon), pc_mean(:,3), pc_min(:,3), pc_max(:,3), 'LineStyle', '--', 'Color', [0.4940 0.1840 0.5560], 'Marker', 'd');
    plot((0:time_horizon), r*ones(time_horizon+1,1),'k-', 'LineWidth', 2);
    drawnow()
    hold off
    xlabel('Time Step, k');
    ylabel('Distance (in meters)');
    axis([-1 time_horizon+1 5 1.1*axis_max])
    set(gca, 'OuterPosition', [0.01,0.01,0.98,0.45]);
    
    l = legend([p1,p2,p3,p4,p5], {'$||A-B||$','$||A-C||$','$||B-C||$',strcat('R=', num2str(r)), 'DC Lower Bound'},...
        'interpreter', 'latex',...
        'Orientation','horizontal');
    set(l,'Position', [0.05,0.94,0.9,0.04],'Units', 'normalized');
    
    set(fig.Children, ...
        'FontName',     'Times', ...
        'FontSize',     20);
    
end