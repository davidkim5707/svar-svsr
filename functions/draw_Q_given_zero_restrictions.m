function [Qdraw] = draw_Q_given_zero_restrictions(Q_init, W, ZF_for_Q, w)
% DRAW_Q_GIVEN_ZERO_RESTRICTIONS
% Construct columns of Q sequentially to satisfy zero/sign restrictions.
%
% Inputs
%   Q_init   : ny x ny matrix (typically zeros(ny) on entry; previous cols may be prefilled)
%   W        : 1 x ny cell, W{j} is (#sign/linear rows for shock j) x ny
%   ZF_for_Q : 1 x ny cell, ZF_for_Q{j} is (#zero-IRF rows for shock j) x ny
%   w        : ny^2 x 1 vector of random weights (e.g., vec(eye(ny)) or randn(ny^2,1))
%
% Outputs
%   Q        : ny x ny orthonormal matrix built column-by-column
%   Qdraw    : Q' (transpose), sometimes convenient for downstream code
%
% Notes
%   - Follows the Rubio-Ramírez–Waggoner–Zha style draw using QR on the
%     stacked restriction matrix's transpose.
%   - For each shock j, we stack:
%         [ Q(:,1:j-1)';  ZF_for_Q{j};  W{j} ]
%     take a QR on its transpose, fix signs on R's diagonal, and select the
%     last s columns where s = size(W{j},1), then combine with random weights.

    % --- Basic checks
    if ~iscell(W) || ~iscell(ZF_for_Q)
        error('W and ZF_for_Q must be 1xny cell arrays.');
    end

    [ny, ny2] = size(Q_init);
    if ny ~= ny2
        error('Q_init must be square (ny x ny).');
    end

    if numel(W) ~= ny || numel(ZF_for_Q) ~= ny
        error('W and ZF_for_Q must have one cell per shock (1 x ny).');
    end

    % if numel(w) ~= ny*ny
    %     error('w must be a vector of length ny^2.');
    % end

    Q = Q_init;
    k = 0;

    for j = 1:ny
        % Rows contributed by W{j} determine how many columns we'll select post-QR
        s = size(W{j}, 1);

        % Pull the random weights for this block
        wj = w(k + (1:s));

        % Stack restrictions for shock j
        %   - previous Q columns enforce orthogonality
        %   - ZF_for_Q{j} encodes IRF-based zero restrictions
        %   - W{j} encodes sign/other linear restrictions (row-wise)
        Mj_tilde = [Q(:, 1:max(j-1,0))';  ZF_for_Q{j};  W{j}];

        % QR on the transpose to get orthonormal basis K
        [K, R] = qr(Mj_tilde');

        % Ensure positive diagonal of R on the block we'll use
        % (this stabilizes the orientation of K)
        for ii = ny - s + 1 : ny
            if R(ii, ii) < 0
                K(:, ii) = -K(:, ii);
            end
        end

        % Take the last s columns as the admissible subspace for shock j
        Kj     = K(:, ny - s + 1 : ny);

        % Form column j of Q by mixing with weights
        Q(:, j) = Kj * wj;

        % Advance the weight offset
        k = k + s;
    end

    % Return both Q and its transpose, since some code expects Q'
    Qdraw = Q;
end
