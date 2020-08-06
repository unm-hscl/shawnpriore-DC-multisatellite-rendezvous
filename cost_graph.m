%% Cost over recent itterations
% Plot costs
figure();
c = min(k-1, 30);
hold on
plot((k-c:k), input_cost(k-c+1:k+1),'k-');
plot((k-c:k), lambda_sum(k-c+1:k+1),'b-');
plot((k-c:k), total_cost(k-c+1:k+1),'r-');
drawnow()
hold off

% Fancy labels
title('Costs Over Recent Itterations', 'interpreter', 'latex');
xlabel('Itteration, $k$', 'interpreter', 'latex');
ylabel('Cost', 'interpreter', 'latex');
legend({'Input Cost','Slack Cost','Total Cost'}, 'interpreter', 'latex');

% Center graph
m = 1.1*max(0.1, max(total_cost(2:end)));
axis([k-c-1 k+1 0 m])