%% Ex2_small_gyro_clean.m
% Example 2 of the paper:  M = I,  K = kappa_K*I,  G = [0 -kappa_G; kappa_G 0],
% kappa_K = kappa_G = 1e-3,  delta M = delta K = 0,  |dG_ij| <= eps |G_ij|.
%
% Produces
%   Table 3  : B_tan,i (22) and the standard tan(Theta) bound (40) for x_1, x_2
%              at eps = 1e-7, with the corresponding gaps; a separate random
%              perturbation dG is drawn for each eigenvector
%   Figure 2 : bound ratio_2 = standard bound / B_tan,i  for
%              eps = 1e-4 ... 1e-9
%   accuracy : residuals, gaps, ||r||/gap
%
% Since delta I = 0 here, (A_L~, I) is a standard Hermitian problem and
% Theorem 4.3 (standard tan Theta bound) applies to the Hermitian matrix A_L.
%
% !! TWO VERSIONS of the standard bound are computed -- see std_variant below:
%   'thm43'  : Theorem 4.3 applied literally to A_L:
%                R = A_L*x~ - x~*(x~'*A_L*x~),  delta = dist(x~'*A_L*x~, {mu_j}_{j~=i})
%              (gap in mu = 1/lambda, about 1 in this example)
%   'paper'  : as in the current Table 3 / Figure 2:
%                ||dA_L||_F / min_j |lambda_i - lambda~_j|
%              (gap in lambda, about 0.001 in this example)
%   With 'paper' the bound ratio is exactly  gap_mu/(2*gap_lambda) = 500
%   (the published value 707.1 corresponds to a factor sqrt(2) instead of the
%   factor 2 in (22)).  With 'thm43' the standard bound is of the same order
%   as B_tan.  Please decide which comparison the paper should report.
%
% Output files: Ex2_table3.csv, Ex2_accuracy.csv, Ex2_bound_ratio.mat

clear; clc; close all;
rng(1);

%% ---------------- system -----------------------------------------------
n = 2;  kappa_K = 1e-3;  kappa_G = 1e-3;
M = eye(n);  K = kappa_K*eye(n);  G = [0 -kappa_G; kappa_G 0];

eps_table = 1e-7;                           % Table 3
eps_list  = [1e-4 1e-5 1e-6 1e-7 1e-8 1e-9];% Figure 2
std_variant = 'paper';                      % 'paper' | 'thm43'   (see header)

%% ---------------- unperturbed problem ---------------------------------
LKi = inv(chol(K,'lower'));  LMi = inv(chol(M,'lower'));
AL  = [LKi*(1i*G)*LKi', LKi*M*LMi'; LMi*M*LKi', zeros(n)];  AL = (AL+AL')/2;
[U, D] = eig(AL);  mu = real(diag(D));  [mu, p] = sort(mu,'ascend');  U = U(:,p);
lam = 1./mu;
pos = find(mu > 0);  [~, o] = sort(abs(lam(pos)), 'descend');  idx = pos(o);
% x_1 <-> |lambda| = 0.03212 ... adjust 'descend'/'ascend' to match the paper's labelling
fprintf('Eigenvalues lambda (QEP: +-1i*lambda): %s\n', mat2str(abs(lam(idx)).', 5));
fprintf('mu = 1/lambda                         : %s\n', mat2str(mu(idx).', 6));

res_unp = vecnorm(AL*U - U*diag(mu));
gap_unp = arrayfun(@(i) min(abs(mu(i) - mu([1:i-1 i+1:end]))), 1:2*n);

%% ---------------- Table 3 ------------------------------------------------
[Btan, Bstd, gap_mu, gap_lam, acc] = one_perturbation(G, AL, U, mu, idx, eps_table, LKi, std_variant);
T3 = table((1:2)', abs(lam(idx)), Btan, gap_mu, Bstd, gap_lam, Bstd./Btan, ...
    'VariableNames', {'i','abs_lambda','Btan_22','gap_mu_for_22','std_bound_40','gap_lambda','ratio2'});
disp(['TABLE 3   (standard bound variant: ' std_variant ')']); disp(T3);
writetable(T3, 'Ex2_table3.csv');

TA = table((1:2)', res_unp(idx)', gap_unp(idx)', res_unp(idx)'./gap_unp(idx)', acc.res_pert, acc.gap_pert, acc.res_pert./acc.gap_pert, ...
    'VariableNames', {'i','res_unpert','gap_unpert','angle_err_unpert','res_pert','gap_pert','angle_err_pert'});
disp('ACCURACY OF COMPUTED EIGENPAIRS  (angle error <= ||r||/gap)'); disp(TA);
writetable(TA, 'Ex2_accuracy.csv');

%% ---------------- Figure 2 : ratio vs eps --------------------------------
ratio2 = zeros(numel(eps_list), 2);
for t = 1:numel(eps_list)
    [Bt_, Bs_] = one_perturbation(G, AL, U, mu, idx, eps_list(t), LKi, std_variant);
    ratio2(t,:) = (Bs_ ./ Bt_).';
end
ef_kvoc1_1 = ratio2(:,1).';  ef_kvoc2_1 = ratio2(:,2).';
save('Ex2_bound_ratio.mat', 'eps_list', 'ef_kvoc1_1', 'ef_kvoc2_1', 'std_variant');
disp('bound ratio_2 for eps_list (rows) and x_1, x_2 (columns):'); disp([eps_list.' ratio2]);

figure; hold on; grid on; box on;
plot(eps_list, ratio2(:,1), 'o', 'LineWidth', 2, 'MarkerSize', 5);
plot(eps_list, ratio2(:,2), 's', 'LineWidth', 2, 'MarkerSize', 5);
set(gca, 'XScale', 'log');
xlabel('random perturbation size \epsilon'); ylabel('Bound ratio_2');
legend('bound ratio for x_1', 'bound ratio for x_2', 'Location', 'best');

%% ---------------- save figures (PNG + EPS) -------------------------------
figs = findobj('Type','figure');
for fh = reshape(figs,1,[])
    set(fh,'PaperPositionMode','auto');
    saveas(fh, sprintf('Ex2_fig%d.png', fh.Number));
    saveas(fh, sprintf('Ex2_fig%d.eps', fh.Number), 'epsc');
end

%% ======================= local functions ==================================
function [Btan, Bstd, gap_mu, gap_lam, acc] = one_perturbation(G, AL, U, mu, idx, eps_pert, LKi, std_variant)
n = size(G,1);  N = 2*n;
k = numel(idx);
Btan = zeros(k,1); Bstd = zeros(k,1); gap_mu = zeros(k,1); gap_lam = zeros(k,1);
res_pert = zeros(k,1); gap_pert = zeros(k,1);
lam = 1./mu;
for q = 1:k
    % a separate random perturbation is drawn for each eigenvector x_q
    R  = rand(n);  dG = eps_pert*((R+R.')/2).*G;        % skew-symmetric, |dG_ij| <= eps|G_ij|
    dAL = blkdiag(LKi*(1i*dG)*LKi', zeros(n));  dAL = (dAL+dAL')/2;
    ALt = AL + dAL;                                      % delta I = 0
    [Ut, Dt] = eig(ALt);  mut = real(diag(Dt));  [mut, p] = sort(mut,'ascend');  Ut = Ut(:,p);
    lamt = 1./mut;
    ndALF = norm(dAL,'fro');                             % alpha = 1, ||Ut||_2 = 1 (delta I = 0)
    i = idx(q);  oth = [1:i-1, i+1:N];
    gap_mu(q)  = min(abs(mut(i) - mu(oth)));
    gap_lam(q) = min(abs(lam(i) - lamt(oth)));
    Btan(q)    = sqrt(2)*ndALF/gap_mu(q);                      % (22) with alpha = 1, ||Ut|| = 1, dI = 0
    switch std_variant
        case 'thm43'
            x  = Ut(:,i);  A1 = x'*AL*x;  Rres = AL*x - x*A1;
            Bstd(q) = norm(Rres) / min(abs(A1 - mu(oth)));
        case 'paper'
            Bstd(q) = ndALF / gap_lam(q);
    end
    r = ALt*Ut(:,i) - mut(i)*Ut(:,i);
    res_pert(q) = norm(r);  gap_pert(q) = min(abs(mut(i) - mut(oth)));
end
acc = struct('res_pert', res_pert, 'gap_pert', gap_pert);
end
