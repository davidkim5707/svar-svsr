function ZF_for_Q = build_ZF_rows_for_Q(info0, A0, lambda_in, tempvar, sign_regime_dependent)
% BUILD_ZF_ROWS_FOR_Q
% Rows are e_i' R(h) A0^{-1} Λ_r^{1/2}, one per zero restriction (and per regime if requested).
% info0.details(d) fields:
%   .var, .shock, .tset     (required)
%   .regimes  (optional) vector of regime indices to bind the zero in; if absent and
%              sign_regime_dependent==1, the zero binds in ALL regimes in lambda_in.
%
% Inputs
%   A0   : n x n (contemporaneous matrix)
%   lambda_in : n x R (R>=1) matrix of regime variances (diagonal entries per regime),
%               or n x 1 vector for a single regime, or [] if not used.
%   tempvar   : structure used by your impulsdtrf (must contain VAR dynamics)
%   sign_regime_dependent :
%         0 -> homoskedastic-style rows (no Λ^{1/2} scaling)
%         1 -> heteroskedastic rows stacked across regimes using Λ_r^{1/2}
%
% Output
%   ZF_for_Q : 1 x n cell; ZF_for_Q{j} is (#rows for shock j) x n

    n = size(A0,1);

    % --- Normalize lambda to an n x R matrix (R>=1) ---
    if isempty(lambda_in)
        lambda_mat = ones(n,1);           % fallback (not used when sign_regime_dependent==0)
    else
        lambda_mat = lambda_in;
        if isvector(lambda_mat); lambda_mat = lambda_mat(:); end
        assert(size(lambda_mat,1) == n, 'lambda_in must have n rows.');
    end
    R = size(lambda_mat,2);

    % --- Precompute variance-free IRFs once with A0^{-1} (no Λ scaling) ---
    smat_base = A0 \ eye(n);              % A0^{-1}, numerically stable
    if ~isfield(info0,'details') || isempty(info0.details)
        ZF_for_Q = repmat({zeros(0,n)}, 1, n);
        return
    end
    Hmax = max(arrayfun(@(d) max(d.tset), info0.details));  % 1 -> R(0), 2 -> R(1), ...
    ir3d_base = impulsdtrf(tempvar, smat_base, Hmax);       % n x Hmax x n

    % --- Init output cells (one per shock j) ---
    ZF_for_Q = repmat({zeros(0,n)}, 1, n);

    % --- Build rows ---
    for d = 1:numel(info0.details)
        iVar   = info0.details(d).var;
        jShock = info0.details(d).shock;
        hs     = info0.details(d).tset(:).';

        % basic checks
        assert(all(hs>=1 & hs<=Hmax), 'tset horizons out of range.');
        assert(iVar>=1 && iVar<=n && jShock>=1 && jShock<=n, 'var/shock index out of range.');

        for h = hs
            % variance-free row for horizon h (1->R(0))
            row_base = reshape(ir3d_base(iVar, h, :), 1, n);  % 1 x n

            if R==1 || sign_regime_dependent == 0
                % Homoskedastic-style: no Λ scaling
                ZF_for_Q{jShock} = [ZF_for_Q{jShock}; row_base];
            else
                % Heteroskedastic: stack one row per regime using Λ_r^{1/2}
                for r = 1:R
                    Dhalf_r = diag(sqrt(max(lambda_mat(:,r),0)));   % Λ_r^{1/2}
                    ZF_for_Q{jShock} = [ZF_for_Q{jShock}; row_base * Dhalf_r];
                end
            end
        end
    end
end
