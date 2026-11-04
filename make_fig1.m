%% make_fig1.m -- Figure 1 (a),(b) from Ex1_bound_ratio.mat, style of the old figure1.m
load('Ex1_bound_ratio.mat');            % ratio_x1, ratio_x2, ratio_x3, seed_used
x = 1:numel(ratio_x1);
pad = @(v) [floor(min(v)*20)/20, ceil(max(v)*20)/20];   % round limits to 0.05

% (a) x_1 and x_2
f1 = figure('Units','centimeters','Position',[2 2 11 11]);
hold on; grid on; box on;
plot(x, ratio_x1, 'or', 'MarkerSize', 3);
plot(x, ratio_x2, 'sb', 'MarkerSize', 3);
xlabel('random perturbation'); ylabel('Bound ratio_1');
legend('eigenvector x_1', 'eigenvector x_2', 'Location', 'northeast');
yl = pad([ratio_x1; ratio_x2]); ylim(yl); xlim([0 numel(x)]);
set(gca, 'XTick', 0:100:numel(x), 'YTick', yl(1):0.05:yl(2), 'FontSize', 10);
axis square;
set(f1,'PaperPositionMode','auto');
saveas(f1, 'fig1a_rev.eps', 'epsc');  saveas(f1, 'fig1a_rev.png');

% (b) x_3
f2 = figure('Units','centimeters','Position',[14 2 11 11]);
hold on; grid on; box on;
plot(x, ratio_x3, '^g', 'MarkerSize', 3);
xlabel('random perturbation'); ylabel('Bound ratio_1');
legend('eigenvector x_3', 'Location', 'northeast');
yl = pad(ratio_x3); ylim(yl); xlim([0 numel(x)]);
set(gca, 'XTick', 0:100:numel(x), 'YTick', yl(1):0.05:yl(2), 'FontSize', 10);
axis square;
set(f2,'PaperPositionMode','auto');
saveas(f2, 'fig1b_rev.eps', 'epsc');  saveas(f2, 'fig1b_rev.png');
disp('fig1 done');
