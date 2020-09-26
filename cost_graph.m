%% Cost over recent itterations
m = 1.1*max(0.1, max(total_cost(2:end)));

% Plot costs
fig = figure();
fig.Units    = 'inches';
fig.Position = [0,1,10,9];

hold on
p1=plot((1:k), log10(input_cost(2:k+1)),'k-');
p2=plot((1:k), log10(lambda_sum(2:k+1)),'b-');
p3=plot((1:k), log10(total_cost(2:k+1)),'r-');
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