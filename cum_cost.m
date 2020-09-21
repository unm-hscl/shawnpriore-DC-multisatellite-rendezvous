norm_holder = zeros(time_horizon, 3);

for i = 1:time_horizon
    norm_holder(i,1) = norm(U_a((i-1)*2 + (1:2)))^2;
    norm_holder(i,2) = norm(U_b((i-1)*2 + (1:2)))^2;
    norm_holder(i,3) = norm(U_c((i-1)*2 + (1:2)))^2;
end

figure()
hold on

plot((0:9),cumsum(norm_holder(:,1)),'k-o');
plot((0:9),cumsum(norm_holder(:,2)),'r-o');
plot((0:9),cumsum(norm_holder(:,3)),'b-o');
drawnow()
hold off

xlabel('Time step, $k$', 'interpreter', 'latex');
ylabel('Cumulative cost', 'interpreter', 'latex');
legend({'$||U_A(t)||^2_2$','$||U_B(t)||^2_2$','$||U_C(t)||^2_2$'},...
    'interpreter', 'latex',...
    'Location','northwest');