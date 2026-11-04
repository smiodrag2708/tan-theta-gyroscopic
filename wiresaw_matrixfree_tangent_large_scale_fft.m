%% wiresaw_matrixfree_tangent_large_scale.m  (FFT version)
%
% CHANGES vs. previous version:
%   * G*x via FFT, O(N log N): G = 2v (Toeplitz + Hankel) diag(1..N)
%   * structured relative perturbation  G~ = S G S,  S = I + eps*diag(d),
%     |dG_jk| <= (2eps+eps^2)|G_jk|, G~ skew-symmetric, also FFT-applicable
%   * ||L_K^{-1} dG L_K^{-*}||_F computed once, blockwise (no N-by-N storage)
%
% Matrix-free large-scale experiment for the wiresaw gyroscopic QEP.
%
% The code computes ONLY the new tangent perturbation bound.
%
% For each dimension N, it computes:
%
%   - k_interest eigenpairs of interest,
%   - m_buffer additional guard eigenvalues,
%   - perturbation norms,
%   - ||Upsilon_tilde||_2 bound without forming the full eigenvector matrix,
%   - complementary spectral gap,
%   - tangent bound B_tan,
%   - Ritz residuals,
%   - separate time for the evaluation of B_tan,
%   - computational times.
%
% At the end, the code produces:
%
%   PaperTable
%   ModesTable
%   DiagnosticsTable
%
% and writes them to CSV files.
%
% -------------------------------------------------------------------------
% Transformed eigenvalue problem:
%
%       A_L u = mu u,             mu = 1/lambda.
%
% Perturbed problem:
%
%       Atilde_L utilde = mutilde Itilde utilde.
%
% Standard transformed perturbed problem:
%
%       D Atilde_L D z = mutilde z,
%
%       D = Itilde^{-1/2},
%       utilde = D z.
%
% IMPORTANT:
%
% eigs(...,'lr') computes the largest positive mu.
% Since mu = 1/lambda, these correspond to the smallest positive lambda,
% i.e. the lowest-frequency vibration modes.
%
% No full N-by-N matrix G and no full 2N-by-2N linearized matrices
% are stored.
%
% =========================================================================

clear;
clc;


%% ========================================================================
%  USER PARAMETERS
% =========================================================================

% Problem dimensions used in the large-scale experiment.
%
% Start with:
N_values = [1000 2500 5000 10000 20000 50000 100000 1000000];
bandWidth = 201;   % half-bandwidth for ||L_K^{-1} dG L_K^{-*}||_F: exact band + rigorous tail, O(N*b)

% If the computation for N = 5000 is sufficiently fast, try:
%
% N_values = [1000 2500 5000 10000];


v        = 0.98;
eps_pert = 1e-8;


% Number of eigenpairs for which the bound is reported
k_interest = 5;


% Additional eigenvalues used as guard eigenvalues for estimating
% the complementary spectral gap
m_buffer = 6;

k_total = k_interest + m_buffer;


% Block size in the matrix-free multiplication by G
blockSize = 500;


% eigs parameters
tol_outer   = 1e-12;
maxit_outer = 300;


%% ========================================================================
%  STORAGE FOR RESULTS
% =========================================================================

numN = length(N_values);


Btan_all = NaN(numN,k_interest);
gap_all  = NaN(numN,k_interest);

lambda_all  = NaN(numN,k_interest);
lambdat_all = NaN(numN,k_interest);

mu_all  = NaN(numN,k_interest);
mut_all = NaN(numN,k_interest);

res_unpert_all = NaN(numN,k_interest);
res_pert_all   = NaN(numN,k_interest);

nearest_index_all = NaN(numN,k_interest);


normdAL_all = NaN(numN,1);
normdI_F_all = NaN(numN,1);
normdI_2_all = NaN(numN,1);

alpha_all = NaN(numN,1);

normUpsilonBound_all = NaN(numN,1);

orthErr_all = NaN(numN,1);


time_norms_all = NaN(numN,1);
time_unpert_all = NaN(numN,1);
time_pert_all = NaN(numN,1);
time_eigs_all = NaN(numN,1);
time_bound_all = NaN(numN,1);
time_total_all = NaN(numN,1);


%% ========================================================================
%  LOOP OVER PROBLEM DIMENSIONS
% =========================================================================

for rr = 1:numN

    N = N_values(rr);

    fprintf('\n\n');
    fprintf('============================================================\n');
    fprintf(' LARGE-SCALE WIRESAW EXPERIMENT\n');
    fprintf(' N = %d, linearization size = %d\n',N,2*N);
    fprintf(' v = %.4f, epsilon = %.1e\n',v,eps_pert);
    fprintf('============================================================\n');


    totalTimer = tic;


    %% ====================================================================
    %  BASIC CHECK
    % =====================================================================

    if k_total >= 2*N
        error('k_total must be smaller than 2*N.');
    end


    %% ====================================================================
    %  DIAGONAL DATA -- O(N) STORAGE
    % =====================================================================

    j_idx = (1:N)';


    % Wiresaw mass matrix:
    %
    % M = 1/2 I

    M_diag = 0.5 * ones(N,1);


    % Wiresaw stiffness matrix:
    %
    % K_jj = j^2*pi^2*(1-v^2)/2

    K_diag = ...
        j_idx.^2 * pi^2 * (1-v^2) / 2;


    % Cholesky factors

    LM_diag = sqrt(M_diag);
    LK_diag = sqrt(K_diag);

    LMinv_diag = 1 ./ LM_diag;
    LKinv_diag = 1 ./ LK_diag;


    % Odd/even index sets needed for G

    isOdd = mod(j_idx,2)==1;

    kOdd  = j_idx(isOdd);
    kEven = j_idx(~isOdd);


    fprintf('Diagonal data constructed -- O(N) storage.\n');


    %% ====================================================================
    %  RANDOM PERTURBATIONS
    % =====================================================================

    % Same random seed for reproducibility

    rng(1);


    dK_diag = ...
        eps_pert * (2*rand(N,1)-1) .* K_diag;


    dM_diag = ...
        eps_pert * (2*rand(N,1)-1) .* M_diag;


    % Structured perturbation of G:  G~ = S G S,  S = I + eps*diag(dG_d)

    dG_d   = 2*rand(N,1) - 1;
    s_cong = 1 + eps_pert*dG_d;

    Gfun  = wiresawG_fft(N,v);                 % y = G*x
    Gtfun = @(x) s_cong .* Gfun(s_cong .* x);  % y = G~*x


    %% ====================================================================
    %  PERTURBED I MATRIX
    % =====================================================================

    % delta I is diagonal in this wiresaw implementation

    dI_diag_full = [ ...
        (LKinv_diag.^2).*dK_diag ; ...
        (LMinv_diag.^2).*dM_diag ];


    Ie_diag = 1 + dI_diag_full;


    if any(Ie_diag <= 0)

        error(['The perturbed matrix Itilde is not positive definite ' ...
               'for N = %d.'],N);

    end


    Dinvsqrt = ...
        1 ./ sqrt(Ie_diag);


    %% ====================================================================
    %  PERTURBATION NORMS
    % =====================================================================

    fprintf('\nComputing perturbation norms...\n');

    timerNorms = tic;


    % Contribution of deltaH = i deltaG to ||delta A_L||_F^2

    S_dG = ...
        dG_frobenius_sq_band(N,v,eps_pert,dG_d,bandWidth);


    % Contribution of deltaM to both off-diagonal blocks of delta A_L

    S_dM_offdiag = ...
        2 * sum( ...
        (LKinv_diag .* dM_diag .* LMinv_diag).^2 );


    normdAL_F = ...
        sqrt(S_dG + S_dM_offdiag);


    % Since delta I is diagonal:

    normdI_F = ...
        norm(dI_diag_full,2);


    normdI_2 = ...
        max(abs(dI_diag_full));


    if normdI_2 >= 0.5

        error(['Assumption ||delta I||_2 < 1/2 is violated ' ...
               'for N = %d.'],N);

    end


    alpha_fac = ...
        sqrt( ...
        (1-normdI_2) / ...
        (1-2*normdI_2) );


    % Bound from the new remark:
    %
    % ||Upsilon_tilde||_2
    %
    % <= 1/sqrt(1-||delta I||_2)
    %
    % Therefore the full diagonalizing matrix is NOT needed.

    normUpsilonBound = ...
        1 / sqrt(1-normdI_2);


    t_norms = toc(timerNorms);


    fprintf('||delta A_L||_F        = %.8e\n',normdAL_F);
    fprintf('||delta I||_F          = %.8e\n',normdI_F);
    fprintf('||delta I||_2          = %.8e\n',normdI_2);
    fprintf('alpha                  = %.12f\n',alpha_fac);
    fprintf('bound for ||Utilde||_2 = %.12f\n',normUpsilonBound);
    fprintf('Norm computation time  = %.2f s\n',t_norms);


    %% ====================================================================
    %  MATRIX-FREE OPERATORS
    % =====================================================================

    matvecAL = ...
        @(x) AL_matvec(x,N,Gfun,LKinv_diag,LM_diag);


    matvecAeL_std = ...
        @(x) AeL_standard_matvec(x,N,Gtfun, ...
        LKinv_diag,LMinv_diag,LM_diag,dM_diag,Dinvsqrt);


    %% ====================================================================
    %  EIGS OPTIONS
    % =====================================================================

    opts.tol    = tol_outer;
    opts.maxit  = maxit_outer;

    opts.isreal = false;
    opts.issym  = true;


    %% ====================================================================
    %  UNPERTURBED EIGENVALUE PROBLEM
    % =====================================================================

    fprintf('\nUnperturbed partial eigensolver...\n');

    timerUnpert = tic;


    [V_part,D_part] = ...
        eigs(matvecAL,2*N,k_total,'lr',opts);


    t_unpert = toc(timerUnpert);


    mu_part = ...
        diag(D_part);


    if max(abs(imag(mu_part))) > 1e-7

        warning(['Unperturbed eigenvalues contain non-negligible ' ...
                 'imaginary parts for N = %d.'],N);

    end


    mu_part = ...
        real(mu_part);


    if any(mu_part <= 0)

        warning('Some computed unperturbed mu values are not positive.');

    end


    lambda_part = ...
        1 ./ mu_part;


    % Sort by increasing positive lambda

    [lambda_part,ord] = ...
        sort(lambda_part,'ascend');


    mu_part = ...
        mu_part(ord);


    V_part = ...
        V_part(:,ord);


    fprintf('Unperturbed eigs time = %.2f s\n',t_unpert);


    %% ====================================================================
    %  PERTURBED EIGENVALUE PROBLEM
    % =====================================================================

    fprintf('Perturbed partial eigensolver...\n');

    timerPert = tic;


    [Vt_part_std,Dt_part] = ...
        eigs(matvecAeL_std,2*N,k_total,'lr',opts);


    t_pert = toc(timerPert);


    mut_part = ...
        diag(Dt_part);


    if max(abs(imag(mut_part))) > 1e-7

        warning(['Perturbed eigenvalues contain non-negligible ' ...
                 'imaginary parts for N = %d.'],N);

    end


    mut_part = ...
        real(mut_part);


    if any(mut_part <= 0)

        warning('Some computed perturbed mu values are not positive.');

    end


    lambdat_part = ...
        1 ./ mut_part;


    [lambdat_part,ordt] = ...
        sort(lambdat_part,'ascend');


    mut_part = ...
        mut_part(ordt);


    Vt_part_std = ...
        Vt_part_std(:,ordt);


    % Return eigenvectors to generalized problem

    Vt_part = ...
        bsxfun(@times,Dinvsqrt,Vt_part_std);


    fprintf('Perturbed eigs time   = %.2f s\n',t_pert);


    %% ====================================================================
    %  ORTHOGONALITY DIAGNOSTIC
    % =====================================================================

    IeVt = ...
        bsxfun(@times,Ie_diag,Vt_part);


    IeOrthCheck = ...
        Vt_part' * IeVt;


    orth_err = ...
        max(max(abs( ...
        IeOrthCheck-eye(k_total) )));


    fprintf('Itilde-orthogonality error = %.4e\n',orth_err);


    %% ====================================================================
    %  RITZ RESIDUALS
    % =====================================================================

    % These residuals provide an a posteriori measure of the accuracy
    % of the partial eigenvalue computation.

    res_unpert = ...
        zeros(k_total,1);


    res_pert = ...
        zeros(k_total,1);


    for jj = 1:k_total

        uj = ...
            V_part(:,jj);


        rj = ...
            matvecAL(uj) ...
            - mu_part(jj)*uj;


        res_unpert(jj) = ...
            norm(rj,2);


        utj = ...
            Vt_part_std(:,jj);


        rtj = ...
            matvecAeL_std(utj) ...
            - mut_part(jj)*utj;


        res_pert(jj) = ...
            norm(rtj,2);

    end


    fprintf('max unperturbed Ritz residual = %.4e\n', ...
            max(res_unpert));

    fprintf('max perturbed Ritz residual   = %.4e\n', ...
            max(res_pert));


    %% ====================================================================
    %  NEW TANGENT BOUND
    % =====================================================================
    %
    % The timer below measures ONLY the numerical evaluation of the
    % tangent bound for the k_interest modes. Printing and diagnostics
    % are deliberately performed after toc, so they are not included
    % in t_bound.
    % =====================================================================

    Btan = ...
        zeros(k_interest,1);

    gapTan = ...
        zeros(k_interest,1);

    nearestIndex = ...
        zeros(k_interest,1);


    timerBound = tic;


    for ii = 1:k_interest

        % Perturbed reciprocal eigenvalue

        mu_ti = ...
            mut_part(ii);


        % Complementary unperturbed eigenvalues among the computed
        % spectral window

        other_idx = ...
            setdiff(1:k_total,ii);


        % Spectral separation:
        %
        % gap_i =
        %
        % min_{j ~= i}
        % | 1/lambda_tilde_i - 1/lambda_j |
        %
        % = min |mu_tilde_i - mu_j|

        gap_values = ...
            abs(mu_ti - mu_part(other_idx));


        [gapTan(ii),locgap] = ...
            min(gap_values);


        nearestIndex(ii) = ...
            other_idx(locgap);


        % ||Lambda_tilde_1||_2
        %
        % = |1/lambda_tilde_i|
        % = |mu_tilde_i|

        normLambdaTilde1 = ...
            abs(mu_ti);


        % New tangent bound:
        %
        % B_tan =
        %
        % 2 alpha ||Upsilon_tilde||_2
        %
        % * ( ||delta A_L||_F
        %     + ||delta I||_F |1/lambda_tilde_i| )
        %
        % / gap_i

        numerator_tan = ...
            2 * ...
            alpha_fac * ...
            normUpsilonBound * ...
            ( ...
              normdAL_F ...
              + normdI_F*normLambdaTilde1 ...
            );


        Btan(ii) = ...
            numerator_tan / gapTan(ii);

    end


    t_bound = toc(timerBound);


    %% ====================================================================
    %  TANGENT-BOUND OUTPUT AND GUARD DIAGNOSTICS
    % =====================================================================

    fprintf('\n');
    fprintf('------------------------------------------------------------\n');
    fprintf(' TANGENT BOUND\n');
    fprintf('------------------------------------------------------------\n');

    fprintf('%-3s %-14s %-14s %-14s %-14s\n', ...
        'i', ...
        'lambda_i', ...
        'lambda_tilde', ...
        'gap_i', ...
        'B_tan');


    for ii = 1:k_interest

        fprintf('%-3d %-14.6e %-14.6e %-14.6e %-14.6e\n', ...
            ii, ...
            lambda_part(ii), ...
            lambdat_part(ii), ...
            gapTan(ii), ...
            Btan(ii));


        % Guard-eigenvalue diagnostic

        if nearestIndex(ii) == k_total

            fprintf(['WARNING for mode %d: nearest complementary ' ...
                     'eigenvalue is the last guard eigenvalue.\n'],ii);

            fprintf('Increase m_buffer and repeat the computation.\n');

        end

    end


    fprintf('Tangent-bound evaluation time = %.6e s\n',t_bound);


    %% ====================================================================
    %  STORE RESULTS
    % =====================================================================

    normdAL_all(rr) = ...
        normdAL_F;


    normdI_F_all(rr) = ...
        normdI_F;


    normdI_2_all(rr) = ...
        normdI_2;


    alpha_all(rr) = ...
        alpha_fac;


    normUpsilonBound_all(rr) = ...
        normUpsilonBound;


    orthErr_all(rr) = ...
        orth_err;


    time_norms_all(rr) = ...
        t_norms;


    time_unpert_all(rr) = ...
        t_unpert;


    time_pert_all(rr) = ...
        t_pert;


    time_eigs_all(rr) = ...
        t_unpert + t_pert;


    time_bound_all(rr) = ...
        t_bound;


    for ii = 1:k_interest

        Btan_all(rr,ii) = ...
            Btan(ii);


        gap_all(rr,ii) = ...
            gapTan(ii);


        lambda_all(rr,ii) = ...
            lambda_part(ii);


        lambdat_all(rr,ii) = ...
            lambdat_part(ii);


        mu_all(rr,ii) = ...
            mu_part(ii);


        mut_all(rr,ii) = ...
            mut_part(ii);


        res_unpert_all(rr,ii) = ...
            res_unpert(ii);


        res_pert_all(rr,ii) = ...
            res_pert(ii);


        nearest_index_all(rr,ii) = ...
            nearestIndex(ii);

    end


    time_total_all(rr) = ...
        toc(totalTimer);


    fprintf('\nTotal time for N = %d: %.2f s\n', ...
            N,time_total_all(rr));


end


%% ========================================================================
%  FINAL TABLE 1: COMPACT TABLE FOR THE PAPER
% =========================================================================
%
% This is the table that is probably most useful in the manuscript.
%
% It reports the first mode, together with computational times.
%
% EigsTime allows us to show that the partial eigenvalue computation
% is the dominant computational cost.
% =========================================================================

PaperTable = table( ...
    N_values(:), ...
    2*N_values(:), ...
    Btan_all(:,1), ...
    gap_all(:,1), ...
    time_norms_all, ...
    time_bound_all, ...
    time_eigs_all, ...
    time_total_all, ...
    'VariableNames', ...
    { ...
    'n', ...
    'LinearizedSize', ...
    'Btan1', ...
    'Gap1', ...
    'NormTime_s', ...
    'BoundTime_s', ...
    'EigsTime_s', ...
    'TotalTime_s' ...
    } ...
    );


fprintf('\n\n');
fprintf('============================================================\n');
fprintf(' TABLE FOR THE MANUSCRIPT\n');
fprintf('============================================================\n');

disp(PaperTable);


%% ========================================================================
%  FINAL TABLE 2: ALL FIVE MODES
% =========================================================================

numberRows = ...
    numN*k_interest;


n_column = ...
    zeros(numberRows,1);

mode_column = ...
    zeros(numberRows,1);

lambda_column = ...
    zeros(numberRows,1);

lambdat_column = ...
    zeros(numberRows,1);

gap_column = ...
    zeros(numberRows,1);

Btan_column = ...
    zeros(numberRows,1);

resU_column = ...
    zeros(numberRows,1);

resT_column = ...
    zeros(numberRows,1);

nearest_column = ...
    zeros(numberRows,1);


row = 0;


for rr = 1:numN

    for ii = 1:k_interest

        row = row + 1;


        n_column(row) = ...
            N_values(rr);


        mode_column(row) = ...
            ii;


        lambda_column(row) = ...
            lambda_all(rr,ii);


        lambdat_column(row) = ...
            lambdat_all(rr,ii);


        gap_column(row) = ...
            gap_all(rr,ii);


        Btan_column(row) = ...
            Btan_all(rr,ii);


        resU_column(row) = ...
            res_unpert_all(rr,ii);


        resT_column(row) = ...
            res_pert_all(rr,ii);


        nearest_column(row) = ...
            nearest_index_all(rr,ii);

    end

end


ModesTable = table( ...
    n_column, ...
    mode_column, ...
    lambda_column, ...
    lambdat_column, ...
    gap_column, ...
    Btan_column, ...
    resU_column, ...
    resT_column, ...
    nearest_column, ...
    'VariableNames', ...
    { ...
    'n', ...
    'Mode', ...
    'Lambda', ...
    'LambdaTilde', ...
    'Gap', ...
    'Btan', ...
    'Residual', ...
    'ResidualTilde', ...
    'NearestGuardIndex' ...
    } ...
    );


fprintf('\n');
fprintf('============================================================\n');
fprintf(' DETAILED RESULTS FOR THE FIRST FIVE MODES\n');
fprintf('============================================================\n');

disp(ModesTable);


%% ========================================================================
%  FINAL TABLE 3: DIAGNOSTICS
% =========================================================================

DiagnosticsTable = table( ...
    N_values(:), ...
    normdAL_all, ...
    normdI_F_all, ...
    normdI_2_all, ...
    alpha_all, ...
    normUpsilonBound_all, ...
    orthErr_all, ...
    time_unpert_all, ...
    time_pert_all, ...
    time_bound_all, ...
    time_total_all, ...
    'VariableNames', ...
    { ...
    'n', ...
    'NormDeltaAL_F', ...
    'NormDeltaI_F', ...
    'NormDeltaI_2', ...
    'Alpha', ...
    'UpsilonNormBound', ...
    'OrthogonalityError', ...
    'UnpertEigsTime_s', ...
    'PertEigsTime_s', ...
    'BoundTime_s', ...
    'TotalTime_s' ...
    } ...
    );


fprintf('\n');
fprintf('============================================================\n');
fprintf(' DIAGNOSTIC TABLE\n');
fprintf('============================================================\n');

disp(DiagnosticsTable);


%% ========================================================================
%  SAVE RESULTS TO CSV FILES
% =========================================================================

writetable( ...
    PaperTable, ...
    'wiresaw_tangent_paper_table.csv');


writetable( ...
    ModesTable, ...
    'wiresaw_tangent_modes_table.csv');


writetable( ...
    DiagnosticsTable, ...
    'wiresaw_tangent_diagnostics.csv');


fprintf('\nResults saved as:\n');

fprintf('  wiresaw_tangent_paper_table.csv\n');
fprintf('  wiresaw_tangent_modes_table.csv\n');
fprintf('  wiresaw_tangent_diagnostics.csv\n');


fprintf('\n============================================================\n');
fprintf(' COMPUTATION FINISHED\n');
fprintf('============================================================\n');


%% ========================================================================
%  LOCAL FUNCTIONS
%  ALL FUNCTIONS MUST REMAIN AT THE END OF THE SCRIPT
% =========================================================================

function Gfun = wiresawG_fft(N,v)
% y = G*x in O(N log N):
%   G_jk = 4v jk/(j^2-k^2) = 2v k [ 1/(j-k) + 1/(j+k) ],  j+k odd,
%   i.e. G = 2v (T + H) diag(1,...,N), T Toeplitz, H Hankel.
L  = 2^nextpow2(4*N);
d  = (-(N-1):(N-1))';  t = zeros(2*N-1,1); o = mod(d,2)~=0; t(o) = 1./d(o);
s  = (2:2*N)';         h = zeros(2*N-1,1); o = mod(s,2)~=0; h(o) = 1./s(o);
Tt = fft(t,L);  Th = fft(h,L);  kk = (1:N)';
Gfun = @(x) Gfft_apply(x,Tt,Th,kk,N,v,L);
end

function y = Gfft_apply(x,Tt,Th,kk,N,v,L)
w = kk.*x;
a = ifft(Tt.*fft(w,L));
b = ifft(Th.*fft(flipud(w),L));
y = 2*v*(a(N:2*N-1) + b(N:2*N-1));
if isreal(x), y = real(y); end
end

%% ========================================================================

function y = AL_matvec(x,N,Gfun,LKinv_diag,LM_diag)
x1 = x(1:N);  x2 = x(N+1:end);
Hz = 1i*Gfun(LKinv_diag.*x1);                       % H = iG
y  = [LKinv_diag.*Hz + LKinv_diag.*LM_diag.*x2 ; ...
      LM_diag.*LKinv_diag.*x1];
end

%% ========================================================================

function y = AeL_matvec(x,N,Gtfun,LKinv_diag,LMinv_diag,LM_diag,dM_diag)
x1  = x(1:N);  x2 = x(N+1:end);
Hz  = 1i*Gtfun(LKinv_diag.*x1);                     % H~ = iG~
off = LKinv_diag.*LM_diag + LKinv_diag.*dM_diag.*LMinv_diag;
y   = [LKinv_diag.*Hz + off.*x2 ; off.*x1];
end

%% ========================================================================

function y = AeL_standard_matvec(x,N,Gtfun, ...
    LKinv_diag,LMinv_diag,LM_diag,dM_diag,Dinvsqrt)
% D Atilde_L D,  D = Itilde^{-1/2}
y = Dinvsqrt .* AeL_matvec(Dinvsqrt.*x,N,Gtfun, ...
    LKinv_diag,LMinv_diag,LM_diag,dM_diag);
end

%% ========================================================================

function s2 = dG_frobenius_sq_band(N,v,eps_pert,dd,b)
% UPPER BOUND for ||L_K^{-1} dG L_K^{-*}||_F^2, dG = S G S - G, in O(N*b):
%   exact sum over |j-k| <= b  +  rigorous tail bound for |j-k| > b:
%   |e_jk| <= eta = 2eps+eps^2,  sum_{|j-k|>b} 1/(j^2-k^2)^2 <= 1/(2 b^2).
% Hence B_tan computed with this value remains an upper bound.
c   = 8*v/(pi^2*(1-v^2));
eta = 2*eps_pert + eps_pert^2;
s2  = 0;
for d = 1:2:min(b,N-1)                 % j-k = d odd  (<=> j+k odd)
    k  = (1:N-d)';  j = k + d;
    e  = eps_pert*(dd(j)+dd(k)) + eps_pert^2*dd(j).*dd(k);
    s2 = s2 + 2*sum( e.^2 ./ (j.^2-k.^2).^2 );   % factor 2: (j,k) and (k,j)
end
s2 = c^2*s2 + c^2*eta^2/(2*b^2);
end

%% ========================================================================

function s2 = dG_frobenius_sq_cong(N,v,eps_pert,dd,blockSize)
% ||L_K^{-1} dG L_K^{-*}||_F^2  for  dG = S G S - G,  S = I + eps*diag(dd).
% Entry (j+k odd):  8v/(pi^2(1-v^2)) * e_jk/(j^2-k^2),
%   e_jk = eps(d_j+d_k) + eps^2 d_j d_k   (no cancellation).
% O(N^2) flops ONCE, O(N*blockSize) memory.
c  = 8*v/(pi^2*(1-v^2));
k  = 1:N;  s2 = 0;
for jstart = 1:blockSize:N
    jb = (jstart:min(jstart+blockSize-1,N))';
    E  = eps_pert*bsxfun(@plus,dd(jb),dd(:)') + eps_pert^2*(dd(jb)*dd(:)');
    D2 = bsxfun(@minus,jb.^2,k.^2).^2;
    W  = (mod(bsxfun(@plus,jb,k),2)==1) ./ D2;      % 0/0 on diagonal
    W(~isfinite(W)) = 0;
    s2 = s2 + c^2*sum(sum(W.*E.^2));
end
end
