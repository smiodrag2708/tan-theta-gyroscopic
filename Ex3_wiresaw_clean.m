%% Ex3_wiresaw_clean.m
% Example 3 of the paper: wiresaw1 problem (NLEVP), moderate dimensions.
%
% Produces
%   Figure 3 : bound ratio_1 = B_sin,i (32) / B_tan,i (22) for the eigenvectors
%              x_1..x_5 of the five smallest positive eigenvalues, v = 0.8,
%              eps = 1e-8, n = 100, 200, 400, 600, 1000
%   Table 4  : lambda_1, B_sin,1 and B_tan,1 for v = 0.3, 0.5, 0.7, 0.9, n = 500
%   accuracy : residuals, gaps and ||r||/gap of the computed eigenpairs
%
% Method (same as Ex1_particle_clean.m):
%   Hermitian transformed linearization (A_L, I), perturbed pair (A_L~, I~)
%   with the unperturbed Cholesky factors; eigenpairs by eig / eig(.,.,'chol');
%   perturbations dX = eps*(R+R')/2 .* X, R ~ U[0,1], so |dX_ij| <= eps|X_ij|.
%
% wiresaw1:  M = I/2,  K = diag(j^2 pi^2 (1-v^2)/2),
%            G_jk = 4 v j k/(j^2-k^2) for j+k odd, 0 otherwise.
%
% Output files: Ex3_fig3_ratio.mat, Ex3_table4.csv, Ex3_accuracy.csv

clear; clc; close all;
seed_used = 852684860;  rng(seed_used);   % seed used for the results in the paper; use rng('shuffle') for a new random perturbation
fprintf('rng seed for this run: %d   (repeat this run with rng(%d))\n', seed_used, seed_used);

eps_pert    = 1e-8;
RG2_variant = 'paper';          % 'paper' (Theorem 4.1) | 'old' (as in tan_theta_ziro_wiresaw.m)
k_vec       = 5;                % eigenvectors of the k smallest positive eigenvalues

%% ---------------- Figure 3 -------------------------------------------------
v_fig   = 0.95;
n_list  = [100 200 400 600 1000];
ratio   = zeros(numel(n_list), k_vec);
acc_all = [];
for t = 1:numel(n_list)
    n = n_list(t);
    [M, G, K] = wiresaw_matrices(n, v_fig);
    [Btan, Bsin, ~, lam1, acc] = wiresaw_bounds(M, G, K, eps_pert, k_vec, RG2_variant);
    ratio(t,:) = (Bsin ./ Btan).';
    acc_all = [acc_all; repmat(n, k_vec, 1), (1:k_vec)', acc.lam, acc.res_unp, acc.gap_unp, acc.res_unp./acc.gap_unp, acc.res_pert, acc.gap_pert, acc.res_pert./acc.gap_pert]; %#ok<AGROW>
    fprintf('n = %5d  lambda_1 = %.6f   ratio_1(x_1..x_5) = %s\n', n, lam1, mat2str(ratio(t,:), 4));
end
save('Ex3_fig3_ratio.mat', 'n_list', 'ratio', 'v_fig', 'eps_pert', 'seed_used');

TA = array2table(acc_all, 'VariableNames', {'n','i','lambda','res_unpert','gap_unpert','angle_err_unpert','res_pert','gap_pert','angle_err_pert'});
disp('ACCURACY OF COMPUTED EIGENPAIRS  (angle error <= ||r||/gap)'); disp(TA);
writetable(TA, 'Ex3_accuracy.csv');

figure; hold on; grid on; box on;
mk = {'o','s','d','^','.'};
for t = 1:numel(n_list)
    plot(1:k_vec, ratio(t,:), mk{t}, 'LineWidth', 1, 'MarkerSize', 8 + 7*(t==5));
end
set(gca, 'YScale', 'log', 'XTick', 1:k_vec, 'XTickLabel', {'x_1','x_2','x_3','x_4','x_5'});
yline(1, 'k:');
xlabel('eigenvectors'); ylabel('Bound ratio_1');
legend(arrayfun(@(n) sprintf('n=%d', n), n_list, 'UniformOutput', false), 'Location', 'best');
xlim([0.5 k_vec+0.5]);

%% ---------------- Table 4 --------------------------------------------------
n_tab  = 500;
v_list = [0.3 0.5 0.7 0.9 0.99];
T4 = zeros(3, numel(v_list));
for t = 1:numel(v_list)
    [M, G, K] = wiresaw_matrices(n_tab, v_list(t));
    [Btan, Bsin, ~, lam1] = wiresaw_bounds(M, G, K, eps_pert, 1, RG2_variant);
    T4(:,t) = [lam1; Bsin(1); Btan(1)];
end
T4tab = array2table(T4, 'VariableNames', arrayfun(@(v) sprintf('v_%g', v), v_list, 'UniformOutput', false), ...
                        'RowNames', {'lambda_1', 'Bsin_1 (32)', 'Btan_1 (22)'});
disp('TABLE 4'); disp(T4tab);
writetable(T4tab, 'Ex3_table4.csv', 'WriteRowNames', true);

%% ---------------- save figures (PNG + EPS) -------------------------------
figs = findobj('Type','figure');
for fh = reshape(figs,1,[])
    set(fh,'PaperPositionMode','auto');
    saveas(fh, sprintf('Ex3_fig%d.png', fh.Number));
    saveas(fh, sprintf('Ex3_fig%d.eps', fh.Number), 'epsc');
end

%% ======================= local functions ==================================
function [M, G, K] = wiresaw_matrices(n, v)
% identical to nlevp('wiresaw1',n,v)
j = (1:n)';  k = 1:n;
G = 8*v*(j*k) ./ (j.^2 - k.^2 + eye(n));      % +eye avoids 0/0 on the diagonal
G(mod(j+k, 2) == 0) = 0;                      % only j+k odd
M = eye(n)/2;  G = G/2;  K = diag(j.^2*pi^2*(1-v^2))/2;
end

function dX = rel_pert(X, eps_pert)
R  = rand(size(X));
dX = eps_pert * ((R + R.')/2) .* X;           % symmetric for M,K; skew-symmetric for G
end

function [Btan, Bsin, s_ex, lam1, acc] = wiresaw_bounds(M, G, K, eps_pert, k_vec, RG2_variant)
n = size(M,1);  N = 2*n;
LKi = inv(chol(K,'lower'));  LMi = inv(chol(M,'lower'));
AL  = [LKi*(1i*G)*LKi', LKi*M*LMi'; LMi*M*LKi', zeros(n)];  AL = (AL+AL')/2;
[U, D] = eig(AL);  mu = real(diag(D));  [mu, p] = sort(mu,'ascend');  U = U(:,p);
lam = 1./mu;

% perturbation
dM = rel_pert(M, eps_pert);  dK = rel_pert(K, eps_pert);  dG = rel_pert(G, eps_pert);
ALt = [LKi*(1i*(G+dG))*LKi', LKi*(M+dM)*LMi'; LMi*(M+dM)*LKi', zeros(n)];  ALt = (ALt+ALt')/2;
It  = blkdiag(eye(n) + LKi*dK*LKi', eye(n) + LMi*dM*LMi');  It = (It+It')/2;
dAL = ALt - AL;  dI = It - eye(N);
[Ut, Dt] = eig(ALt, It, 'chol');  mut = real(diag(Dt));  [mut, p] = sort(mut,'ascend');  Ut = Ut(:,p);
lamt = 1./mut;

% eigenvectors of the k smallest positive eigenvalues lambda  (= largest positive mu)
pos = find(mu > 0);  [~, o] = sort(lam(pos), 'ascend');  idx = pos(o(1:k_vec));
lam1 = lam(idx(1));

nI2 = norm(dI);  alpha = sqrt((1-nI2)/(1-2*nI2));
nUt = norm(Ut);  kUt = cond(Ut);
ndALF = norm(dAL,'fro');  ndIF = norm(dI,'fro');  nAinvdA = norm(AL\dAL,'fro');

Btan = zeros(k_vec,1); Bsin = zeros(k_vec,1); s_ex = zeros(k_vec,1);
res_unp = zeros(k_vec,1); gap_unp = zeros(k_vec,1); res_pert = zeros(k_vec,1); gap_pert = zeros(k_vec,1);
for q = 1:k_vec
    i = idx(q);  oth = [1:i-1, i+1:N];
    gap_mu  = min(abs(mut(i) - mu(oth)));
    Btan(q) = (1/cos(pi/4))*alpha*nUt*(ndALF + ndIF*abs(mut(i)))/gap_mu;
    RG1 = min(abs(lam(i) - lamt(oth)) ./ abs(lam(oth)));
    switch RG2_variant
        case 'paper', RG2 = min(abs(lam(i) - lamt(oth))) / abs(lamt(i));
        case 'old',   RG2 = min(abs(lam(i) - lamt(oth)) ./ abs(lam(oth)));
    end
    Bsin(q) = kUt*(nAinvdA/RG2 + ndIF/RG1);
    s_ex(q) = norm(U(:,oth)'*Ut(:,i));
    res_unp(q)  = norm(AL*U(:,i) - mu(i)*U(:,i));
    gap_unp(q)  = min(abs(mu(i) - mu(oth)));
    res_pert(q) = norm(ALt*Ut(:,i) - mut(i)*(It*Ut(:,i))) / norm(Ut(:,i));
    gap_pert(q) = min(abs(mut(i) - mut(oth)));
end
acc = struct('lam', lam(idx), 'res_unp', res_unp, 'gap_unp', gap_unp, 'res_pert', res_pert, 'gap_pert', gap_pert);
end
