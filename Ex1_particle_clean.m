%% Ex1_particle_clean.m
% Example 1 of the paper: particle on a rotating ring (Datta et al.), n = 3.
%
% Produces
%   Table 1 : B_tan,i (22), B_sin,i (32) and the exact value of |sin angle|
%             for x_1, x_2, x_3 (lambda = 6.0860i, 3.1895i, 0.8878i)
%   Table 2 / Figure 1 : bound ratio_1 = B_sin,i / B_tan,i for 500 random
%             perturbations
%   accuracy check : Ritz residuals, gaps and ||r||/gap for the computed
%             eigenpairs (paragraph "Accuracy of the computed eigenvectors")
%
% Method (as described in "Computational setup"):
%   * Hermitian transformed linearization  A_L = L^{-1} A L^{-*},  I,
%     perturbed pair  (A_L~, I~),  I~ = I + delta I,  with the SAME Cholesky
%     factors L_M, L_K of the unperturbed M, K.
%   * eigenvalues mu = 1/lambda and eigenvectors by eig:
%       [U ,D ] = eig(AL)             -> U unitary
%       [Ut,Dt] = eig(ALt,It,'chol')  -> Ut  It-orthonormal (Ut'*It*Ut = I)
%   * perturbations  dX = eps*(R+R')/2 .* X,  R ~ U[0,1],  so that
%       |dX_ij| <= eps |X_ij|   (condition (37) with eps),
%     dM, dK symmetric, dG skew-symmetric.
%
% Output files:  Ex1_table1.csv, Ex1_table2.csv, Ex1_accuracy.csv,
%                Ex1_bound_ratio.mat (ratio_x1, ratio_x2, ratio_x3)

clear; clc; close all;
% Fixed seed reproduces the tables and figures of the paper.
seed_used = 852683514;  rng(seed_used);   % seed used for the results in the paper; use rng('shuffle') for a new random perturbation
fprintf('rng seed for this run: %d   (repeat this run with rng(%d))\n', seed_used, seed_used);

%% ---------------- system matrices ------------------------------------
w  = [1; 2; 1];  wn = 3;  gam = 1/2;
M  = eye(3);
G  = [    0    -2*w(3)   2*w(2);
       2*w(3)     0     -2*w(1);
      -2*w(2)   2*w(1)      0   ];
K  = [2*wn^2-w(2)^2-w(3)^2,      w(1)*w(2),                 w(1)*w(3);
      w(1)*w(2),                 2*wn^2*gam-w(1)^2-w(3)^2,  w(2)*w(3);
      w(1)*w(3),                 w(2)*w(3),                 2*wn^2*gam-w(1)^2-w(2)^2];
n = size(M,1);

eps_pert = 1e-8;      % perturbation level in (37)
n_rand   = 500;       % number of random perturbations for Figure 1 / Table 2
RG2_variant = 'paper';   % 'paper': RG2 = min |lambda_i-lambda~_j| / |lambda~_i|  (Theorem 4.1)
                         % 'old'  : RG2 = min |lambda_i-lambda~_j| / |lambda_j|    (as in tan_theta_ziro.m)

%% ---------------- unperturbed transformed pencil ---------------------
[AL, LKi, LMi] = transformed_pencil(M, G, K);
[U, mu] = herm_eig(AL);               % mu ascending, U unitary
lam = 1./mu;                          % QEP eigenvalues are  +-1i*|lam|
pos = find(mu > 0);                   % three eigenpairs with positive mu
% paper's ordering: x_1,x_2,x_3 <-> |lambda| = 6.0860, 3.1895, 0.8878  (mu ascending)
[~, ord] = sort(mu(pos), 'ascend');  idx = pos(ord);
fprintf('Eigenvalues lambda (QEP: +-1i*lambda):  %s\n', mat2str(abs(lam(idx)).', 5));

% residuals / gaps of the unperturbed problem
res_unp = vecnorm(AL*U - U*diag(mu));
gap_unp = arrayfun(@(i) min(abs(mu(i) - mu([1:i-1 i+1:end]))), 1:2*n);

%% ---------------- Table 1 : one random perturbation -------------------
[Btan, Bsin, s_ex, t_ex, acc] = one_perturbation(M, G, K, AL, LKi, LMi, U, mu, idx, eps_pert, RG2_variant);

T1 = table((1:3)', abs(lam(idx)), Btan, Bsin, s_ex, t_ex, ...
    'VariableNames', {'i','abs_lambda','Btan_22','Bsin_32','exact_sin','exact_tan'});
disp('TABLE 1'); disp(T1);
writetable(T1, 'Ex1_table1.csv');

TA = table((1:3)', res_unp(idx)', gap_unp(idx)', res_unp(idx)'./gap_unp(idx)', acc.res_pert, acc.gap_pert, acc.res_pert./acc.gap_pert, ...
    'VariableNames', {'i','res_unpert','gap_unpert','angle_err_unpert','res_pert','gap_pert','angle_err_pert'});
disp('ACCURACY OF COMPUTED EIGENPAIRS  (angle error <= ||r||/gap)'); disp(TA);
writetable(TA, 'Ex1_accuracy.csv');

%% ---------------- Figure 1 / Table 2 : 500 random perturbations ------
ratio = zeros(n_rand, 3);
for t = 1:n_rand
    [Bt_, Bs_] = one_perturbation(M, G, K, AL, LKi, LMi, U, mu, idx, eps_pert, RG2_variant);
    ratio(t,:) = (Bs_ ./ Bt_).';           % bound ratio_1 = B_sin / B_tan
end
ratio_x1 = ratio(:,1); ratio_x2 = ratio(:,2); ratio_x3 = ratio(:,3);
save('Ex1_bound_ratio.mat', 'ratio_x1', 'ratio_x2', 'ratio_x3', 'eps_pert', 'seed_used');

T2 = table((1:3)', min(ratio)', max(ratio)', 'VariableNames', {'i','ratio_min','ratio_max'});
disp('TABLE 2  (min / max of bound ratio_1 over 500 perturbations)'); disp(T2);
writetable(T2, 'Ex1_table2.csv');

figure;
subplot(1,2,1); hold on; grid on; box on;
plot(1:n_rand, ratio_x1, 'or', 'MarkerSize', 3);
plot(1:n_rand, ratio_x2, 'sb', 'MarkerSize', 3);
xlabel('random perturbation'); ylabel('Bound ratio_1');
legend('eigenvector x_1', 'eigenvector x_2', 'Location', 'best');
subplot(1,2,2); grid on; box on;
plot(1:n_rand, ratio_x3, '^g', 'MarkerSize', 3);
xlabel('random perturbation'); ylabel('Bound ratio_1');
legend('eigenvector x_3', 'Location', 'best');

%% ---------------- save figures (PNG + EPS) -------------------------------
figs = findobj('Type','figure');
for fh = reshape(figs,1,[])
    set(fh,'PaperPositionMode','auto');
    saveas(fh, sprintf('Ex1_fig%d.png', fh.Number));
    saveas(fh, sprintf('Ex1_fig%d.eps', fh.Number), 'epsc');
end

%% ======================= local functions ==============================
function [AL, LKi, LMi] = transformed_pencil(M, G, K)
% A_L = L^{-1} A L^{-*}  with  L = blkdiag(L_K, L_M),  H = 1i*G  (Hermitian)
n   = size(M,1);
LKi = inv(chol(K, 'lower'));
LMi = inv(chol(M, 'lower'));
AL  = [LKi*(1i*G)*LKi',  LKi*M*LMi';
       LMi*M*LKi',       zeros(n)];
AL  = (AL + AL')/2;                    % enforce exact Hermitian symmetry
end

function [U, mu] = herm_eig(A)
[U, D] = eig(A);  mu = real(diag(D));
[mu, p] = sort(mu, 'ascend');  U = U(:,p);
end

function dX = rel_pert(X, eps_pert)
% |dX_ij| <= eps |X_ij|;  symmetric if X symmetric, skew-symmetric if X skew
R  = rand(size(X));
dX = eps_pert * ((R + R.')/2) .* X;
end

function [Btan, Bsin, s_ex, t_ex, acc] = one_perturbation(M, G, K, AL, LKi, LMi, U, mu, idx, eps_pert, RG2_variant)
n   = size(M,1);  N = 2*n;
dM  = rel_pert(M, eps_pert);  dK = rel_pert(K, eps_pert);  dG = rel_pert(G, eps_pert);
Mt  = M + dM;  Gt = G + dG;
% perturbed pair with the unperturbed Cholesky factors
ALt = [LKi*(1i*Gt)*LKi',  LKi*Mt*LMi';
       LMi*Mt*LKi',       zeros(n)];
ALt = (ALt + ALt')/2;
It  = blkdiag(eye(n) + LKi*dK*LKi',  eye(n) + LMi*dM*LMi');
It  = (It + It')/2;
dAL = ALt - AL;  dI = It - eye(N);

[Ut, Dt] = eig(ALt, It, 'chol');  mut = real(diag(Dt));
[mut, p] = sort(mut, 'ascend');  Ut = Ut(:,p);       % It-orthonormal
lam  = 1./mu;  lamt = 1./mut;

nI2   = norm(dI);  alpha = sqrt((1-nI2)/(1-2*nI2));
nUt   = norm(Ut);  kUt = cond(Ut);
ndALF = norm(dAL, 'fro');  ndIF = norm(dI, 'fro');
nAinvdA = norm(AL\dAL, 'fro');                        % ||A_L^{-1} dA_L||_F  for (32)

k = numel(idx);
Btan = zeros(k,1); Bsin = zeros(k,1); s_ex = zeros(k,1); t_ex = zeros(k,1);
res_pert = zeros(k,1); gap_pert = zeros(k,1);
for q = 1:k
    i = idx(q);  oth = [1:i-1, i+1:N];
    % ---- new bound (22) --------------------------------------------
    gap_mu  = min(abs(mut(i) - mu(oth)));             % |1/lambda~_i - 1/lambda_j|
    Btan(q) = (1/cos(pi/4))*alpha*nUt*(ndALF + ndIF*abs(mut(i))) / gap_mu;
    % ---- sine bound (32), Theorem 4.1 -------------------------------
    RG1 = min(abs(lam(i) - lamt(oth)) ./ abs(lam(oth)));
    switch RG2_variant
        case 'paper', RG2 = min(abs(lam(i) - lamt(oth))) / abs(lamt(i));
        case 'old',   RG2 = min(abs(lam(i) - lamt(oth)) ./ abs(lam(oth)));
    end
    Bsin(q) = kUt*(nAinvdA/RG2 + ndIF/RG1);
    % ---- exact values ------------------------------------------------
    S = U(:,oth)' * Ut(:,i);                          % Ut(:,i) is It-normalized
    s_ex(q) = norm(S);
    t_ex(q) = s_ex(q) / sqrt(max(1 - s_ex(q)^2, 0));
    % ---- accuracy of the perturbed eigenpair -------------------------
    r = ALt*Ut(:,i) - mut(i)*(It*Ut(:,i));
    res_pert(q) = norm(r) / norm(Ut(:,i));
    gap_pert(q) = min(abs(mut(i) - mut(oth)));
end
acc = struct('res_pert', res_pert, 'gap_pert', gap_pert);
end
