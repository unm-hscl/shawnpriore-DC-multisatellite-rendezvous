%% Plot of Motion Path
function motion_path_graph(all_a, all_b, all_c, target_set_a, target_set_b, target_set_c, have_legend)
    % make graph
    figure();
    hold on

    % plot path lines
    for i = 1:size(all_a,2)
        plot(all_a(1:4:end, i), all_a(2:4:end, i), 'k-o');
    end
    for i = 1:size(all_b,2)
        plot(all_b(1:4:end, i), all_b(2:4:end, i), 'b-o');
    end
    for i = 1:size(all_c,2)
        plot(all_c(1:4:end, i), all_c(2:4:end, i), 'r-o');
    end

    % mark starting points and plot target regions
    plot(all_a(1, 1), all_a(2, 1), 'ko', 'MarkerFaceColor', 'k');
    plot( polyshape(Polyhedron(target_set_a.A([1;2;5;6],1:2), target_set_a.b([1;2;5;6])).V),...
        'FaceColor', 'k', ...
        'FaceAlpha',0.1); 

    plot(all_b(1, 1), all_b(2, 1), 'bo', 'MarkerFaceColor', 'b');
    plot( polyshape(Polyhedron(target_set_b.A([1;2;5;6],1:2), target_set_b.b([1;2;5;6])).V),...
        'FaceColor', 'b', ...
        'FaceAlpha',0.1); 
    
    plot(all_c(1, 1), all_c(2, 1), 'ro', 'MarkerFaceColor', 'r');
    plot( polyshape(Polyhedron(target_set_c.A([1;2;5;6],1:2), target_set_c.b([1;2;5;6])).V),...
        'FaceColor', 'r', ...
        'FaceAlpha',0.1); 

    % draw graph
    drawnow()
    hold off

    % Fancy labels
    xlabel('$x$ (in meters)', 'interpreter', 'latex')
    ylabel('$y$ (in meters)', 'interpreter', 'latex')
    if have_legend == 1
        legend({'Vehicle A', 'Vehicle B', 'Vehicle C', 'Initial Location', 'Target Set' }, ...
            'interpreter', 'latex', ...
            'Location', 'northeast');
    end    
end