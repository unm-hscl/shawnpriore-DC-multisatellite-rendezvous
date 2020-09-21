%% Cost over recent itterations
% Plot costs
figure();
hold on
plot((1:k), log10(input_cost(2:k+1)),'k-');
plot((1:k), log10(lambda_sum(2:k+1)),'b-');
plot((1:k), log10(total_cost(2:k+1)),'r-');
drawnow()
hold off

% Fancy labels
xlabel('Itteration, $k$', 'interpreter', 'latex');
ylabel('$\log_{10}(Cost)$', 'interpreter', 'latex');
legend({'Input Cost','Slack Cost','Total Cost'}, 'interpreter', 'latex');

% Center graph
m = 1.1*max(0.1, max(total_cost(2:end)));
axis([0 k+1 -inf inf])