%% Graph of distance at each time point
% Add initial point to mean
all_a = [x_0_a; mean_X_a];
all_b = [x_0_b; mean_X_b];
all_c = [x_0_c; mean_X_c];

% Create empty vector
dist_ab = zeros(time_horizon,1);
dist_ac = zeros(time_horizon,1);
dist_bc = zeros(time_horizon,1);

% Find distances between vehicles using L_\infty norm
for i = 0:time_horizon
   dist_ab(i+1) = norm(all_a((1:2)+(i*4))- all_b((1:2)+(i*4)), Inf);
   dist_ac(i+1) = norm(all_a((1:2)+(i*4))- all_c((1:2)+(i*4)), Inf);
   dist_bc(i+1) = norm(all_b((1:2)+(i*4))- all_c((1:2)+(i*4)), Inf);
end

% Plot graph
figure();
hold on
plot((0:time_horizon), dist_ab,'r-o');
plot((0:time_horizon), dist_ac,'g-o');
plot((0:time_horizon), dist_bc,'b-o');
plot((0:time_horizon), r*ones(time_horizon+1,1),'k-');
drawnow()
hold off

% Fancy labels
title('$L_\infty$ Distance Between Vehicles', 'interpreter', 'latex');
xlabel('Time Step, $t$', 'interpreter', 'latex');
ylabel('Distance (in meters)', 'interpreter', 'latex');
legend({'$||A-B||_\infty$','$||A-C||_\infty$','$||B-C||_\infty$',strcat('R=', num2str(r))}, 'interpreter', 'latex');

% Center graph
axis([-1 time_horizon+1 0 1.1*max([dist_ab;dist_ac;dist_bc])])