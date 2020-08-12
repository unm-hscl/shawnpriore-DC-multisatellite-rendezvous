%% Plot of Motion Path
% Add initial point to mean
all_a = [x_0_a; mean_X_a];
all_b = [x_0_b; mean_X_b];
all_c = [x_0_c; mean_X_c];

% Plot graph
figure();
hold on
plot(all_a(1:4:end), all_a(2:4:end), 'g-o');
plot(all_b(1:4:end), all_b(2:4:end), 'b-o');
plot(all_c(1:4:end), all_c(2:4:end), 'r-o');
drawnow()
hold off

% Fancy labels
title('Motion Path', 'interpreter', 'latex')
xlabel('$x$ (in meters)', 'interpreter', 'latex')
ylabel('$y$ (in meters)', 'interpreter', 'latex')
legend({'Vehicle A','Vehicle B','Vehicle C'}, 'interpreter', 'latex');

% Center graph
m = 150;
axis([-m m -m m])