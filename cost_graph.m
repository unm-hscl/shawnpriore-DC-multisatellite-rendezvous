function cost_graph(input_cost, lambda_sum, total_cost)
    %% Cost over recent itterations
    % Plot costs
    fig = figure();
    fig.Units    = 'inches';
    fig.Position = [0,0,10,9];

    k = size(input_cost,1);
    
    hold on
    p1=plot((1:k), log10(input_cost),'k-');
    p2=plot((1:k), log10(lambda_sum),'b-');
    p3=plot((1:k), log10(total_cost),'r-');
    drawnow()
    hold off

    xlabel('Itteration, k');
    ylabel('$\log_{10}(Cost)$', 'interpreter', 'latex');
    axis([0 k+1 -inf inf])
    set(gca, 'OuterPosition', [0.01, 0.07, 0.98, 0.93]);


    l = legend([p1,p2,p3],{'Input Cost','Slack Cost','Total Cost'},...
        'Orientation','horizontal');
    set(l,'Position', [0.23,0.01,0.54,0.05],'Units', 'normalized');

    set(fig.Children, ...
        'FontName',     'Times', ...
        'FontSize',     20);
end