% figure 
% plot3(Beta_t,Beta_t,ti_MECNOMA21_1, '--o')
% grid on
% xlabel('\beta^t')
% ylabel('\beta^e')
% zlabel('Total execution time (s)')
% 
% figure
% plot3(Beta_t,Beta_t,en_MECNOMA21_1, '--o')
% grid on
% xlabel('\beta^t')
% ylabel('\beta^e')
% zlabel('Total execution energy (J)')

% Create figure
figure('OuterPosition',[61 318 835 482]);

% Create axes
axes1 = axes('Position',...
    [0.137823198749905 0.161282164014256 0.839047615692728 0.79898771403988]);
hold(axes1,'on');

plot(Beta_t,ti_MECNOMA21_1, '-o', 'LineWidth', 2.0, 'MarkerSize', 12, 'Color', [0.00,0.00,1.00] )
grid on
xlabel('\beta^t')
ylabel('Total execution time (s)')
yyaxis right
plot(Beta_t,en_MECNOMA21_1, '-*', 'LineWidth',2.0, 'MarkerSize', 12, 'Color', [0.00,1.00,0.00] )
ylabel('Total execution energy (J)')
legend("time", "energy")

box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
set(axes1,'FontSize',17.5,'PlotBoxAspectRatio',[5 2 1]);