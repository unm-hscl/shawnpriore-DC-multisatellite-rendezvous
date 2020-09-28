%% Plot of Motion Path
function trajectories(all_a, all_b, all_c, all_a_pc, all_b_pc, all_c_pc, target_set_a, target_set_b, target_set_c, particles)
    
    red = [0.4660 0.6740 0.1880];
    blue = [0.3010 0.7450 0.9330];
    green = [0.6350 0.0780 0.1840];
    
    % make graph
    fig = figure();
    fig.Units    = 'inches';
    fig.Position = [0.75,-1,10,11.25];
    
    subplot(2,1,1);
    hold on
    p1 = plot(all_a(1:4:end), all_a(2:4:end), 'Color', red, 'Marker', 'h', 'LineWidth', 1);
    p2 = plot(all_b(1:4:end), all_b(2:4:end), 'Color', blue, 'Marker', 'p', 'LineWidth', 1);
    p3 = plot(all_c(1:4:end), all_c(2:4:end), 'Color', green, 'Marker', '^', 'LineWidth', 1);
    p4 = plot(all_a(1, 1), all_a(2, 1), 'Color', red, 'Marker', 'h', 'MarkerFaceColor', 'k');
    plot(all_b(1, 1), all_b(2, 1), 'Color', blue, 'Marker', 'p', 'MarkerFaceColor', 'k');
    plot(all_c(1, 1), all_c(2, 1), 'Color', green, 'Marker', '^', 'MarkerFaceColor', 'k');
    p5 = plot( polyshape(Polyhedron(target_set_a.A([1;2;5;6],1:2), target_set_a.b([1;2;5;6])).V),...
        'FaceColor', red, ...
        'FaceAlpha',0.1); 
    plot( polyshape(Polyhedron(target_set_b.A([1;2;5;6],1:2), target_set_b.b([1;2;5;6])).V),...
        'FaceColor', blue, ...
        'FaceAlpha',0.1); 
    plot( polyshape(Polyhedron(target_set_c.A([1;2;5;6],1:2), target_set_c.b([1;2;5;6])).V),...
        'FaceColor', green, ...
        'FaceAlpha',0.1); 
    xlabel('x (in meters)')
    ylabel('y (in meters)')
    drawnow()
    axis([-20 120 -15 25])
    hold off
    set(gca, 'OuterPosition', [0.01,0.47,0.98,0.45]);


    subplot(2,1,2);
    hold on
    for i = 1:particles
        plot(all_a_pc(1:4:end,i), all_a_pc(2:4:end,i), 'LineStyle', '--', 'Color', red, 'Marker', 'h', 'LineWidth', 1);
        plot(all_b_pc(1:4:end,i), all_b_pc(2:4:end,i), 'LineStyle', '--', 'Color', blue, 'Marker', 'p', 'LineWidth', 1);
        plot(all_c_pc(1:4:end,i), all_c_pc(2:4:end,i), 'LineStyle', '--', 'Color', green, 'Marker', '^', 'LineWidth', 1);
    end
    plot(all_a_pc(1, 1), all_a_pc(2, 1), 'Color', red, 'Marker', 'h', 'MarkerFaceColor', 'k');
    plot(all_b_pc(1, 1), all_b_pc(2, 1), 'Color', blue, 'Marker', 'p', 'MarkerFaceColor', 'k');
    plot(all_c_pc(1, 1), all_c_pc(2, 1), 'Color', green, 'Marker', '^', 'MarkerFaceColor', 'k');
    plot( polyshape(Polyhedron(target_set_a.A([1;2;5;6],1:2), target_set_a.b([1;2;5;6])).V),...
        'FaceColor', red, ...
        'FaceAlpha',0.1); 
    plot( polyshape(Polyhedron(target_set_b.A([1;2;5;6],1:2), target_set_b.b([1;2;5;6])).V),...
        'FaceColor', blue, ...
        'FaceAlpha',0.1); 
    plot( polyshape(Polyhedron(target_set_c.A([1;2;5;6],1:2), target_set_c.b([1;2;5;6])).V),...
        'FaceColor', green, ...
        'FaceAlpha',0.1); 
    xlabel('x (in meters)')
    ylabel('y (in meters)')
    drawnow()
    axis([-20 120 -15 25])
    hold off
    set(gca, 'OuterPosition', [0.01,0.01,0.98,0.45]);

    l = legend([p1,p2,p3,p4,p5], {'Vehicle A', 'Vehicle B', 'Vehicle C', 'Initial Location', 'Target Set' },...
        'Orientation','horizontal');
    set(l,'Position', [0.05,0.94,0.9,0.04],'Units', 'normalized');

    set(fig.Children, ...
        'FontName',     'Times', ...
        'FontSize',     20);
end