function cum_cost(U_a, U_b, U_c, U_a_pc, U_b_pc, U_c_pc, time_horizon)
	
    red = [0.4660 0.6740 0.1880];
    blue = [0.3010 0.7450 0.9330];
    green = [0.6350 0.0780 0.1840];
    
    norm_holder = zeros(time_horizon, 3);
	norm_holder_pc = zeros(time_horizon, 3);

	for i = 1:time_horizon
	    norm_holder(i,1) = norm(U_a((i-1)*2 + (1:2)))^2;
	    norm_holder(i,2) = norm(U_b((i-1)*2 + (1:2)))^2;
	    norm_holder(i,3) = norm(U_c((i-1)*2 + (1:2)))^2;
	    norm_holder_pc(i,1) = norm(U_a_pc((i-1)*2 + (1:2)))^2;
	    norm_holder_pc(i,2) = norm(U_b_pc((i-1)*2 + (1:2)))^2;
	    norm_holder_pc(i,3) = norm(U_c_pc((i-1)*2 + (1:2)))^2;
    end
    
    max_u = max([sum(norm_holder, 1), sum(norm_holder_pc, 1)]) + 1;
    
	fig = figure();
    fig.Units    = 'inches';
    fig.Position = [0,1,10,9];
    
	hold on
	p1 = plot((0:(time_horizon-1)),cumsum(norm_holder(:,1)), 'Color', red, 'Marker', 'h', 'LineWidth', 1);
	p2 = plot((0:(time_horizon-1)),cumsum(norm_holder(:,2)), 'Color', blue, 'Marker', 'p', 'LineWidth', 1);
	p3 = plot((0:(time_horizon-1)),cumsum(norm_holder(:,3)), 'Color', green, 'Marker', '^', 'LineWidth', 1);	
	plot((0:(time_horizon-1)),cumsum(norm_holder_pc(:,1)), 'LineStyle', '--', 'Color', red, 'Marker', 'h', 'LineWidth', 1);
	plot((0:(time_horizon-1)),cumsum(norm_holder_pc(:,2)), 'LineStyle', '--', 'Color', blue, 'Marker', 'p', 'LineWidth', 1);
	plot((0:(time_horizon-1)),cumsum(norm_holder_pc(:,3)), 'LineStyle', '--', 'Color', green, 'Marker', '^', 'LineWidth', 1);
	drawnow()
	xlabel('Time step, k');
	ylabel('Cumulative cost');
    axis([0 time_horizon-1 0 max_u])
	hold off
    set(gca, 'OuterPosition', [0.01, 0.01, 0.98, 0.93]);


	l = legend([p1,p2,p3], {'$||U_A(k)||^2_2$','$||U_B(k)||^2_2$','$||U_C(k)||^2_2$'},...
	    'interpreter', 'latex',...
        'Orientation','horizontal');
    set(l,'Position', [0.37,0.94,0.26,0.05],'Units', 'normalized');

    set(fig.Children, ...
        'FontName',     'Times', ...
        'FontSize',     20);
end