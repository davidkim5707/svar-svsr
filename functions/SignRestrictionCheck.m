
function [SignCheck, accepted_shock, irf_save, Aplus, Qdraw, x1_rotated, stability, lambda, By, fsign] = SignRestrictionCheck(x1, var1, SignRestrictions, NarrativeRestrictions, lags, lcA0, constant, ny, horizon, sign_regime_dependent, penalty_offdiagonal_on, lcLmd, fix_first, Q_given)
%--------------------------------------------------------------------------
% SignRestrictionCheck
%   Checks whether a draw satisfies user-defined sign restrictions
%
% Inputs:
%   - x1        : current draw of structural parameters (vectorized)
%   - var1      : structure from posterior evaluation (contains Bdraw)
%   - SignRestrictions : cell array of string expressions representing sign restrictions
%   - lags      : number of VAR lags
%   - lcA0      : logical mask for free parameters in A0
%   - constant  : boolean flag (true if VAR includes constant term)
%   - ny        : number of endogenous variables in the VAR

%   - Q         : 1 if using random orthogonal Q; 0 for permutation
%
% Outputs:
%   - SignCheck       : 1 if a shock satisfies the restrictions, 0 otherwise
%   - accepted_shock  : index of the accepted structural shock (always 1 here)
%   - fsign           : the sign pattern of the accepted shock
%   - irf_save        : impulse responses [ny x horizon x 1] for accepted shock
%   - Aplus           : lag coefficient matrices [ny x ny x lags]
%   - Qdraw           : orthogonal matrix used for sign rotation (empty if Q == 0)
%--------------------------------------------------------------------------
    if nargin < 13 || isempty(fix_first)
        fix_first = 0;
    end
    fix_first = logical(fix_first);

    stability      = 1;
    
    if ~isempty(NarrativeRestrictions)
        narrativerestrictions = cell(length(NarrativeRestrictions.time_index), 1);
        for i = 1:length(NarrativeRestrictions.time_index)
            t_idx = NarrativeRestrictions.time_index(i);
            s_idx = NarrativeRestrictions.shock_index(i);
            sgn   = NarrativeRestrictions.sign(i);
            if sgn > 0
                narrativerestrictions{i} = sprintf('v([%d],%d)>0', t_idx, s_idx);
            elseif sgn < 0
                narrativerestrictions{i} = sprintf('v([%d],%d)<0', t_idx, s_idx);
            else
                error('Narrative sign must be +1 or -1');
            end
        end
    else
        narrativerestrictions = {};
    end

    % --- Extract A+ from Bdraw ---
    if constant
        Btrim = var1.Bdraw(1:end-1, :);  % remove constant row
    else
        Btrim = var1.Bdraw;
    end

    Aplus = NaN(ny, ny, lags);
    for iLag = 1:lags
        rows = (iLag - 1) * ny + (1:ny);
        Aplus(:, :, iLag) = Btrim(rows, :)';
    end

    % --- Reconstruct A0 from structural draw ---
    A0 = zeros(ny);
    A0(lcA0) = x1(1:sum(lcA0(:)));

    % --- Extract lambda matrix ---
    if all(lcLmd(:) == 0)  % No regimes case
        lambda  = ones(ny, 1);
        nlambda = 1;
    
    else
        nlambda   = size(lcLmd, 2);              % # regimes
        lambda    = ones(ny, nlambda);           % start with ones (helps when fixing col 1)
    
        % Prevent overwriting col 1 if it's fixed-to-1
        if fix_first
            lcLmd(:, 1) = false;                 % <-- was 1; must be false to keep col1 untouched
            lambda(:, 1) = 1;                    % hard-fix first regime
        else
            lambda(:, :) = 0;                    % if not fixing, we'll fill from params
        end
    
        % Fill only free entries from parameter vector
        offA0 = sum(lcA0(:));
        nFree = sum(lcLmd(:));
        if nFree > 0
            lambda(lcLmd) = x1(offA0 + (1:nFree));
        end
    
        % Last column = residual to enforce row-sum constraint
        if fix_first
            if nlambda >= 3
                lambda(:, end) = (nlambda - 1) - sum(lambda(:, 2:end-1), 2);
            else
                % nlambda == 2 → only one free regime left; it must be 1
                lambda(:, end) = 1;
            end
        else
            lambda(:, end) = nlambda - sum(lambda(:, 1:end-1), 2);
        end
    end

    % --- Prepare VAR dynamics structure ---
    By = zeros(ny, ny, lags);
    for iLag = 1:lags
        By(:, :, iLag) = A0 \ Aplus(:, :, iLag);
    end
    tempvar = By;
    
    % === Stability Check ===
    Phi = reshape(By, ny, ny * lags);  % Companion matrix top block
    if lags > 1
        Companion = [Phi; eye(ny*(lags-1)), zeros(ny*(lags-1), ny)];
    else
        Companion = Phi;
    end

    roots = eig(Companion);
    isStable = all(abs(roots) < 1);

    if ~isStable
        [SignCheck,accepted_shock,irf_save,Aplus,Qdraw,x1_rotated,stability,fsign] = deal(0,[],[],[],[],[],0,[]);
        return
    end

    % --- Identification --- 
    shock_token_cells = cellfun(@(s) regexp(s, ',\s*(\d+)\)', 'tokens', 'once'), ...
        SignRestrictions, 'UniformOutput', false);
    
    shock_ids = cellfun(@(c) sscanf(c{1}, '%d'), shock_token_cells);

    % Get the unique shock numbers
    accepted_shock = unique(shock_ids);
    
    % Check Zero restrictions
    if constant
        zero_info = SetupZeroInfo(ny, 1, SignRestrictions);
    else
        zero_info = SetupZeroInfo(ny, 0, SignRestrictions);
    end
    
    if ~isempty(zero_info)
        % --- Precompute ZF rows for the current regime r ---
        ZF    = build_ZF_rows_for_Q(zero_info, A0, lambda, tempvar, sign_regime_dependent); % 1×nvar cell, rows×nvar
        Q0    = zeros(ny);                 % or carry over partially-filled Q
        if sign_regime_dependent==0
            W = zero_info.W;
        else
            W = build_W_from_QZ(Q0, ZF);      % per draw
        end
        w     = DrawW(W);                        % now consistent
        Qdraw = draw_Q_given_zero_restrictions(Q0, W, ZF, w);
    else
        % --- Random orthogonal matrix Q from Haar measure ---
        if nargin < 14 || isempty(Q_given)
            [Qdraw, R]  = qr(randn(ny));
            In          = diag(sign(diag(R)));
            Qdraw       = Qdraw  * In;
        else
            Qdraw       = Q_given;                   % use provided rotation
        end
    end
    Qdraw = Qdraw';

    % --- Calculate the Q'*A0 for recording the posterior ---
    A0_rotated = Qdraw * A0;    
    
    if nlambda == 1
        % --- No regimes case ---
        smat = A0 \ Qdraw';
        irTemp = impulsdtrf(tempvar, smat, horizon);

        [SignCheck, fsign, ~] = checkrestrictions_per_shock(SignRestrictions, irTemp, 1e-12, true);
        
        % Apply per-shock flips for downstream use/plots:
        irTemp_narrative_regime = irTemp .* reshape(fsign, 1, 1, []);
    else
        NumSignCheck_array = zeros(1, nlambda);
        irTemp_narrative_cell = cell(1, nlambda);
        fsign_cell = cell(1, nlambda);

        for i = 1:nlambda
            lambda_i = lambda(:,i);
            Dhalf = diag(sqrt(max(lambda_i, 0)));
            if sign_regime_dependent == 0
                smat_base = (A0 \ Qdraw') * Dhalf;
            else
                smat_base = (A0 \ Dhalf) * Qdraw';              
            end
        
            irTemp = impulsdtrf(tempvar, smat_base, horizon);
            [SignCheck_regime, fsign_local, ~] = checkrestrictions_per_shock(SignRestrictions, irTemp, 1e-12, true);
            if isfield(NarrativeRestrictions, 'regime_index')
                irTemp_narrative_cell{i} = irTemp .* reshape(fsign_local, 1, 1, []); 
                fsign_cell{i} = fsign_local;  % Track fsign per iteration
            else
                irTemp_narrative_cell{i} = [];
                fsign_cell{i} = fsign_local;
            end
        
            is_penalty = false(length(accepted_shock), 1);
            Qlambda = Qdraw * (diag(lambda_i)) * Qdraw';   % == Q Λ Q'
            for k = 1:length(accepted_shock)
                j = accepted_shock(k);
                own_var = abs(Qlambda(j, j));
                cross_var = abs(Qlambda(j, :));
                cross_var(j) = 0;
                is_penalty(k) = own_var < 3*max(cross_var);
            end
        
            if penalty_offdiagonal_on
                penalty = sum(is_penalty);
            else
                penalty = 0.0;
            end
        
            if SignCheck_regime && penalty < 1
                NumSignCheck_array(i) = 1;
            end
        end
        
        SignCheck = sum(NumSignCheck_array) == nlambda;
        smat = A0 \ Qdraw';
        fsign = fsign_cell{1};
    end    
             
    if ~isempty(narrativerestrictions) && SignCheck
        % Check narrative restrictions using existing function (handles multiple restrictions)
        if nlambda>1
            idx = NarrativeRestrictions.regime_index;  % Use first regime for initial check
            v = var1.udraw_true * Qdraw' * fsign_cell{idx}(1);  % rotate structural shocks
        else
            v = var1.udraw_true * Qdraw' * fsign(1);  % rotate structural shocks
        end
        NarrativeCheck = checkNarrativeRestrictions(narrativerestrictions, v);

        % Process each narrative restriction for dominance testing
        num_restrictions = length(NarrativeRestrictions.time_index);
        dominance_checks_passed = false(num_restrictions, 1);

        for nar_idx = 1:num_restrictions
            % Get regime-specific information for this restriction
            if nlambda > 1
                if isfield(NarrativeRestrictions, 'regime_indices') && length(NarrativeRestrictions.regime_indices) >= nar_idx
                    regime_idx = NarrativeRestrictions.regime_indices(nar_idx);
                else
                    regime_idx = NarrativeRestrictions.regime_index;  % Fallback to single regime
                end
                irTemp_narrative_regime = irTemp_narrative_cell{regime_idx};
                v_regime = var1.udraw_true * Qdraw' * fsign_cell{regime_idx}(1);
            else
                v_regime = var1.udraw_true * Qdraw' * fsign(1);
            end

            % Historical decomposition for this specific restriction
            t_idx = NarrativeRestrictions.time_index(nar_idx);
            s_idx = NarrativeRestrictions.shock_index(nar_idx);

            % Handle variable_index (may be scalar or array)
            if isfield(NarrativeRestrictions, 'variable_index')
                if length(NarrativeRestrictions.variable_index) >= nar_idx
                    v_idx = NarrativeRestrictions.variable_index(nar_idx);
                else
                    v_idx = NarrativeRestrictions.variable_index(1);  % Use first if not enough elements
                end
            else
                v_idx = 1;  % Default to first variable if not specified
            end

            % Perform historical decomposition for this restriction
            HD = getHDs_fast(irTemp_narrative_regime, v_regime(t_idx, :), v_idx);
            abs_contrib_target = abs(HD(s_idx));
            abs_contrib_others = sum(abs(HD)) - abs_contrib_target;

            % Check if target shock is main contributor for this restriction
            dominance_checks_passed(nar_idx) = (abs_contrib_target > abs_contrib_others);
        end

        % Aggregation rule: ALL restrictions must pass both narrative and dominance tests
        is_main_contributor = all(dominance_checks_passed);

        if NarrativeCheck == 1 && is_main_contributor == 1
            SignCheck = 1;
        else
            SignCheck = 0;
        end
    end
    
    % === If sign and (narrative) restrictions are satisfied ===
    if SignCheck
        n_shock = length(accepted_shock);
        irf_save = zeros(ny, horizon, n_shock);
        irTemp = impulsdtrf(tempvar, smat, horizon);
        for j = 1:n_shock
            irf_save(:, :, j) = fsign(j) * irTemp(:, :, accepted_shock(j));
        end
        % === Construct x1_rotated ===
        x1_rotated = x1;
        x1_rotated(1:sum(lcA0(:))) = A0_rotated(lcA0);
    
        if nlambda > 1
            offset_start = sum(lcA0(:));
            flat_lmd_indices = find(lcLmd);
            rotated_lmd_vec = zeros(size(flat_lmd_indices));
    
            count_diag = 0;
            for r = 1:nlambda
                lam_diag = lambda(:, r);
                lam_mat = diag(lam_diag);

                if sign_regime_dependent == 0 % no rotation of structural variances
                    lambda_rot = lam_mat;
                else
                    lambda_rot = Qdraw * lam_mat * Qdraw';
                end
    
                % Save diagonal of rotated Lambda for x1_rotated (except last regime)
                if r < nlambda
                    for i = 1:ny
                        count_diag = count_diag + 1;
                        rotated_lmd_vec(count_diag) = lambda_rot(i, i);
                    end
                end
            end
    
            x1_rotated(offset_start + (1:length(rotated_lmd_vec))) = rotated_lmd_vec;
            x1_rotated = x1_rotated(:);  % Ensure column vector
        end
    else
        irf_save = [];
        x1_rotated = [];
    end

end
