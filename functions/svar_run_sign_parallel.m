function output = svar_run_sign_parallel(data, options)
%==========================================================================
% FUNCTION: svar
%==========================================================================
% PURPOSE:
%   Estimate Bayesian Structural VAR (SVAR) model with optional 
%   regime-switching heteroskedasticity, using Gibbs sampling and 
%   Metropolis-Hastings algorithm.
%
% INPUTS:
%   data    : Struct with fields:
%              - y               : Filtered data (T x nVar)
%              - varnames        : Cell array of variable names
%              - filtered_data   : Full dataset (table with date & vars)
%
%   options : Struct with SVAR settings:
%              - lags, hor, ndraw, nburn, nsep
%              - Minnesota prior hyperparams: minn_prior_*
%              - regimes        : Regime break dates (empty = no regime)
%              - firstobs, heterosked, etc.
%
% OUTPUT:
%   output  : Struct containing posterior draws & estimated parameters
%
%==========================================================================

%% Unpack data & options
y               = data.y;
varnames        = data.varnames;
filtered_data   = data.filtered_data;

ny              = options.ny;
lags            = options.lags;
constant        = options.constant;
ndraw           = options.ndraw;
nburn           = options.nburn;
nsep            = options.nsep;
nobs            = options.nobs;
tvA             = options.tvA;
noLmd           = options.noLmd;
irf_horizon     = options.irf_horizon;

minn_prior_lambda   = options.minn_prior_lambda;
minn_prior_mu       = options.minn_prior_mu;

%********************************************************
%% Estimate mode of A0 and regimes
%********************************************************
% Initial A0 estimate (incorporating A0_restriction)
% Run reduced-form VAR
[T,~] = size(y);
xdata = ones(T,constant);
seedrfmodel  = rfvar3(y, lags, xdata, [], minn_prior_lambda,minn_prior_mu,[],[]);

A0_restriction = options.A0_restriction;
if isempty(A0_restriction)
    lcA0 = tril(true(ny));  % includes diagonal
else
    lcA0 = A0_restriction ~= 0;  % TRUE where nonzero, FALSE where zero
end

% Compute Residual Covariance Matrix
tempAinv = chol((seedrfmodel.u' * seedrfmodel.u) / size(seedrfmodel.u, 1))';
tempA = inv(tempAinv);

% Convert Free Parameters to Vector Form
seedA = tempA;                    
seedA(~lcA0) = 0;  % apply structural zero restrictions

%% Handle regimes
regimes = options.regimes;
nBreaks = length(regimes);          % Number of break dates
breakInd = [];                      % Initialize as empty
Tsigbrk  = 0;

% Check if regimes datetime(2014,10,1)is not empty
if ~isempty(regimes)
    breakInd = NaN(1,nBreaks);     % Preallocate index vector
    
    % Find index of each regime date in data.date
    for i = 1:nBreaks
        idx = find(filtered_data.dates == regimes(i));
        if ~isempty(idx)
            breakInd(i) = idx;
        end
    end
    Tsigbrk = [0,breakInd];         %Adding zero at beginning of Tsigbrk

    % define tparam or alpha k.
    if isfield(options, 'tparam') && ~isempty(options.tparam)     
        tparam = options.tparam;
        tscale = options.tscale;
    else
        alpha = options.alpha;
        K = options.K;
        tparam = [];
        tscale = [];
    end

end

% Construct regime vector (includes start, break indices - lags, and end)
regime_vec = [0, breakInd - lags, nobs];  % Ensure breakInd is row vector
nRegimes = length(regime_vec) - 1;         % Number of regimes

seedLmd = ones(ny, nRegimes);       % Create seedLmd matrix
lcLmd = true(ny, nRegimes);         % Create logical matrix of TRUEs (in MATLAB: true)
lcLmd(:, nRegimes) = false;         % Set the last column to FALSE

% NEW: optionally fix first regime to 1
if isfield(options, 'fix_first_regime') && options.fix_first_regime == 1
    lcLmd(:, 1) = false;
    fix_first_regime = options.fix_first_regime;
else
    fix_first_regime = [];
end

%********************************************************
% Estimate A0 and LMD
%********************************************************
% Prepare for x and Hessian to optimize
seedA_vec = seedA(lcA0);            % Extract elements from seedA where lcA0 is true
seedLmd_vec = seedLmd(lcLmd);       % Extract elements from seedLmd where lcLmd is true

% Construct an x vector, if one is not already provided
if ~tvA  % A0 is constant
    seedx = [seedA_vec; seedLmd_vec];   % Concatenate both vectors: Column vector
else
    if noLmd
        % Estimate without a prior on lmd --- set constant to 1
        seedx = repmat(seedA_vec, nRegimes, 1);
        lcLmd = zeros(ny, 1);         % Create logical matrix of TRUEs (in MATLAB: true)
    else
        % Estimate with the prior on lmd
        seedx = [repmat(seedA_vec, nRegimes, 1); seedLmd_vec];
    end
end

nAparams = sum(lcA0(:));      % Count number of A0 parameters
if tvA
    nAparams = nAparams*nRegimes;
end
nLmdparams = sum(lcLmd(:));   % Count number of Lmd parameters
diag_elements = [ones(nAparams, 1); 1e-5 * ones(nLmdparams, 1)];  % Diagonal elements for seedH
seedH = diag(diag_elements);  % Create diagonal matrix

% Set the parameters for the optimization
crit = 1e-10;
nit  = options.nit;
max_compute =1;

bvar_obj = @(seedx) bvar_posterior(seedx, y, lags, lcA0, lcLmd, Tsigbrk, [], options);

switch max_compute
    case 1 % cminwel
        % First optimization attempt
        [~, xh, ~, H, ~, ~, retcodeh] = csminwel(bvar_obj, seedx, seedH, [], crit, nit, 2);
        term1 = retcodeh;

        % Extra check if termination code is not 1 (not convergence)
        if term1 ~= 1
            [~, xh1, ~, ~, ~, ~, ~] = csminwel(bvar_obj, xh, seedH, [], crit, nit, 2);
            tol = 1e-12;  % Close to 1e4 * .Machine$double.eps
            if any((xh - xh1) < tol)
                xh = xh1;
            end
        end

    case 2 % fmincon with constraints
        optim_options = optimset('display','iter', 'MaxFunEvals',1000000, 'TolFun',1e-4, 'TolX',1e-4, 'Algorithm', 'interior-point');
        [xh, ~, ~, ~, ~, ~, H] = fmincon(bvar_obj, seedx, constraint_A, constraint_b, [], [], [], [], [], optim_options);
end

[~, p] = chol(H); 
if p > 0
    % Not PSD — fix it
    [V, D] = eig((H + H') / 2);       % ensure symmetry
    D(D < 1e-10) = 1e-10;             % set small/negative eigenvalues to small positive
    H = V * D * V';
end

%********************************************************
%% Gibbs sampling and MCMC
%********************************************************
% Set MCMC starting point
% rng('default'); % Set seed to 1234 (or any integer)
x0 = xh + mvnrnd(zeros(length(xh), 1), H, 1)';

Sigma = options.hsnscale*H;
draw_x0 = zeros(ndraw, length(x0));
if constant == 1
    draw_Phi = zeros(ny * lags + 1, ny, ndraw);
else
    draw_Phi = zeros(ny * lags, ny, ndraw);
end
draw_udraw = zeros(T-lags,ny,ndraw);
draw_udraw_true = zeros(T-lags,ny,ndraw);
draw_lh = zeros(ndraw,1);
draw_likelihood = zeros(ndraw,1);
draw_A0prior = zeros(ndraw,1);
draw_lambdaprior = zeros(ndraw,1);

if any(lcLmd(:))  % equivalent to ~iszero(lcLmd)
    draw_dout = zeros(nobs, ny, 1000);
    delta0 = draw_dout(:,:,1);
else
    delta0 = [];
end

% === SIGN RESTRICTION SETUP ===
if isfield(options, 'SignRestrictions') && ~isempty(options.SignRestrictions)

    % Retrieve sign restriction-related options
    SignRestrictions         = options.SignRestrictions;

    if isfield(options, 'NarrativeRestrictions') && ~isempty(options.NarrativeRestrictions)
        NarrativeRestrictions = options.NarrativeRestrictions;
        
        % If time_index is specified, determine regime index for each restriction
        if isfield(NarrativeRestrictions, 'time_index') && ~isempty(NarrativeRestrictions.time_index)
            t_narr_array = NarrativeRestrictions.time_index;

            % Ensure time_index is a row vector for consistent processing
            if iscolumn(t_narr_array)
                t_narr_array = t_narr_array';
            end

            % Initialize regime indices array
            regime_indices = zeros(size(t_narr_array));

            % Process each narrative restriction time index
            for nar_idx = 1:length(t_narr_array)
                t_narr = t_narr_array(nar_idx)+options.lags;

                % Display the date for this restriction
                if t_narr <= length(filtered_data.dates)
                    fprintf('Narrative restriction %d at time %d: %s\n', ...
                        nar_idx, t_narr, datestr(filtered_data.dates(t_narr)));
                end

                % Determine which regime this time index falls into using regime_vec
                regime_found = false;
                for r = 1:nRegimes
                    t_start = regime_vec(r) + 1;
                    t_end   = regime_vec(r+1);
                    if t_narr >= t_start && t_narr <= t_end
                        regime_indices(nar_idx) = r;
                        regime_found = true;
                        fprintf('  -> Assigned to regime %d (periods %d-%d)\n', r, t_start, t_end);
                        break;
                    end
                end

                if ~regime_found
                    warning('Narrative restriction %d at time %d does not fall within any regime', ...
                        nar_idx, t_narr);
                    regime_indices(nar_idx) = 1; % Default to first regime
                end
            end

            % Store regime indices for each narrative restriction
            NarrativeRestrictions.regime_indices = regime_indices;

            % For backward compatibility, also store regime_index as the first restriction's regime
            NarrativeRestrictions.regime_index = regime_indices(1);
        end
    else
        NarrativeRestrictions = [];
    end
    
    sign_regime_dependent    = options.sign_regime_dependent;
    inneriteration           = options.sign_inneriteration;
    penalty_offdiagonal_on   = options.penalty_offdiagonal_on;

    % === Parse shock indices from restriction strings ===
    % Extract the shock number from expressions like 'y(i,1:5,SHOCK)'
    shock_token_cells = cellfun(@(s) regexp(s, ',\s*(\d+)\)', 'tokens', 'once'), ...
                                SignRestrictions, 'UniformOutput', false);
    
    % Convert tokens to integer shock indices
    shock_ids = cellfun(@(c) sscanf(c{1}, '%d'), shock_token_cells);

    % Identify unique shocks affected by sign restrictions
    accepted_shock = unique(shock_ids);

    % === Initialize storage for sign restriction diagnostics ===
    sign_attempts         = zeros(ndraw, 1);        % Attempts before successful draw
    sign_accepted_shock   = zeros(ndraw, 1);        % Accepted shock per draw

    if any(lcLmd(:))
        draw_fsign            = zeros(1, ny, ndraw);        % Final sign flip (+1 or -1)
    end

    % === Initialize storage for IRFs under sign restrictions ===
    % General case: multiple shocks with restrictions
    sign_irfs_temp = zeros(ny, irf_horizon, max(accepted_shock), ndraw);

    % Store orthogonal Q matrices
    sign_Qdraws = zeros(ny, ny, ndraw);
    
    % Store posterior distribution, likelihood, prior density
    draw_lh_wQ = zeros(ndraw,1);
    draw_likelihood_wQ = zeros(ndraw,1);
    draw_A0prior_wQ = zeros(ndraw,1);
    draw_lambdaprior_wQ = zeros(ndraw,1);
    
    % Store Q orthogonal matrix
    draw_Qdraw = zeros(ny, ny, ndraw);
    
    % Store Lambda to get irfs by regime
    draw_Lambda = zeros(ny,nRegimes,ndraw);

    % Store BY
    draw_By = zeros(ny, ny, lags, ndraw);

    % Flag to enable sign restriction filtering
    use_signrestrictions = true;
    
else
    % No restrictions provided
    use_signrestrictions = false;
end


%********************************************************
% Run Gibbs sampling
%********************************************************
%% Burn-in
if nburn > 0

    % if use_signrestrictions
    %     % ---- one-time setup BEFORE the burn-in loop ----
    %     %parpool('Processes');                   % enable getCurrentTask()
    % end

    acceptance_burning = 0;
    burnin_start_time = tic;  % Start global timer
    has_hetero = any(lcLmd(:));

    for iburn = 1:nburn        

        % === CASE 1: Sign restrictions only, no hetero: SKIP M-H ===
        if use_signrestrictions && ~has_hetero
            [lh0, ~, ~, ~, var0] = eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);
            SignCheck = 0;
            while SignCheck == 0
                x1 = x0 + mvnrnd(zeros(length(x0),1), Sigma)';
                [lh1, ~, ~, ~, var1] = eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);

                if isempty(var1) || ~isstruct(var1) || ~isfield(var1, 'Bdraw')
                    continue;
                end

                Ncand = inneriteration;
                Qbatch = draw_Q_batch(ny, Ncand);   % seed varies by iburn (and/or by proposal count)
                pass  = false(Ncand,1);
                stabv = false(Ncand,1);
        
                parfor ii = 1:Ncand
                    Qk = Qbatch(:,:,ii);   % each worker gets a different Q
        
                    [sc, ~, ~, ~, ~, ~, st] = SignRestrictionCheck( ...
                        x1, var1, ...
                        SignRestrictions, NarrativeRestrictions, ...
                        lags, lcA0, constant, ny, irf_horizon, ...
                        sign_regime_dependent, penalty_offdiagonal_on, lcLmd, fix_first_regime, Qk);
        
                    pass(ii)  = (sc == 1);
                    stabv(ii) = (st == 1);
                end
        
                idx = find(pass & stabv, 1, 'first');
                if ~isempty(idx)
                    x0 = x1; var0 = var1; lh0 = lh1;
                    acceptance_burning = acceptance_burning + 1;
                    SignCheck = 1;
                end     
            end

        % === CASE 2: Sign + hetero: use M-H + sign check ===
        elseif use_signrestrictions && has_hetero
            [lh0, ~, ~, ~, var0] = eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
            SignCheck = 0;
            
            while SignCheck == 0
                x1 = x0 + mvnrnd(zeros(length(x0),1), Sigma)';            
                [lh1, ~, ~, ~, var1] = eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
            
                if log(rand()) < (lh0 - lh1)
                    if isempty(var1) || ~isstruct(var1) || ~isfield(var1,'Bdraw')
                        continue;
                    end
            
                    Ncand = inneriteration;
                    Qbatch = draw_Q_batch(ny, Ncand);   % seed varies by iburn (and/or by proposal count)
                    pass  = false(Ncand,1);
                    stabv = false(Ncand,1);
            
                    parfor ii = 1:Ncand
                        Qk = Qbatch(:,:,ii);   % each worker gets a different Q
            
                        [sc, ~, ~, ~, ~, ~, st] = SignRestrictionCheck( ...
                            x1, var1, ...
                            SignRestrictions, NarrativeRestrictions, ...
                            lags, lcA0, constant, ny, irf_horizon, ...
                            sign_regime_dependent, penalty_offdiagonal_on, lcLmd, fix_first_regime, Qk);
            
                        pass(ii)  = (sc == 1);
                        stabv(ii) = (st == 1);
                    end
            
                    idx = find(pass & stabv, 1, 'first');
                    if ~isempty(idx)
                        x0 = x1; var0 = var1; lh0 = lh1;
                        acceptance_burning = acceptance_burning + 1;
                        SignCheck = 1;
                    else
                        %fprintf("300 Qdraws fail to satisfy the restrictions, so draw another x1.\n")
                    end
                end
            end

        % === CASE 3: No sign restrictions: standard M-H ===
        else
            if has_hetero
                [lh0, ~, ~, ~, var0] = eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
            else
                [lh0, ~, ~, ~, var0] = eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);
            end

            x1 = x0 + mvnrnd(zeros(length(x0),1), Sigma)';

            if has_hetero
                [lh1, ~, ~, ~, var1] = eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
            else
                [lh1, ~, ~, ~, var1] = eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);
            end

            if log(rand()) < lh0 - lh1
                x0 = x1; var0 = var1; lh0 = lh1;
                acceptance_burning = acceptance_burning + 1;
            end
        end

        % === Delta update if heteroskedastic ===
        if has_hetero
            eout = var0.udraw(1:size(delta0, 1), :) .* exp(delta0);
            if isempty(tparam)
                delta0 = drawdelta(eout, alpha, K);
            else
                delta0 = drawt(eout, tparam, (tscale^2) * tparam);
            end
        end

        time_to_eachdraw = toc(burnin_start_time);
        fprintf('[INFO] Time to reach %4d accepted burn-ins: %.2f seconds\n', iburn,time_to_eachdraw);
        % fprintf("Accepted burnings are %d\n", acceptance_burning);
        % === Progress report ===
        if mod(iburn, 10) == 0 || iburn == nburn
            fprintf('Burn-in iteration: %4d / %4d | Acceptance: %.2f%%\n', ...
                iburn, nburn, 100 * acceptance_burning / iburn);
            time_to_10 = toc(burnin_start_time);
            fprintf('[INFO] Time to reach %4d accepted burn-ins: %.2f seconds\n', iburn,time_to_10);
        end
    end

    fprintf('Acceptance number in burn-in: %d\n', acceptance_burning);
end

fprintf("Sampling process starts\n");
acceptance_sampling = 0;
samplein_start_time = tic;  % Start global timer

for idraw = 1:ndraw
    for isep = 1:nsep
        if isep==5
            fprintf("inner sampling is %d\n", isep)
        end

        % === STEP 1: Evaluate current posterior ===
        delta_input = [];
        if any(lcLmd(:))
            delta_input = delta0;
        end
        [lh0, likelihood0, A0_prior0, lambda_prior0, var0] = eval_posterior(x0, y, lags, lcA0, lcLmd, Tsigbrk, delta_input, tvA, options);

        % -------------------------------
        % Sign-restricted + Heteroskedasticity case
        % -------------------------------
        if use_signrestrictions && any(lcLmd(:))
        
            SignCheck = 0;
            nDrawAttempts = 0;
    
            while SignCheck == 0
                nDrawAttempts = nDrawAttempts + 1;
                x1 = x0 + mvnrnd(zeros(length(x0),1), Sigma)';
                
                [lh1, likelihood1, A0_prior1, lambda_prior1, var1] = eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);

                if log(rand()) < lh0 - lh1

                    if isempty(var1) || ~isstruct(var1) || ~isfield(var1, 'Bdraw')
                        fprintf("Valid var1 is not detected\n");
                        continue;
                    end

                    % --- pre-draw distinct Q's on the client (vary seed by iburn and/or proposal count)
                    Ncand  = inneriteration;
                    Qbatch = draw_Q_batch(ny, Ncand);   % your helper with optional seed
                    
                    % --- preallocate outputs that vary by ii
                    pass      = false(Ncand,1);
                    stabv     = false(Ncand,1);
                    irf_list  = cell(Ncand,1);
                    Qdraw_list = cell(Ncand,1);
                    xrot_list = cell(Ncand,1);
                    lambda_list = cell(Ncand,1);
                    By_list = cell(Ncand, 1);
                    fsign_list = cell(Ncand, 1);

                    parfor ii = 1:Ncand
                        Qk = Qbatch(:,:,ii);   % deterministic, different per ii
                    
                        [sc, ~, irf_tmp, ~, Qdraw_tmp, xrot_tmp, st, lambda_tmp, By_tmp, fsign_tmp] = SignRestrictionCheck( ...
                            x1, var1, ...
                            SignRestrictions, NarrativeRestrictions, ...
                            lags, lcA0, constant, ny, irf_horizon, ...
                            sign_regime_dependent, penalty_offdiagonal_on, lcLmd, fix_first_regime, Qk);
                    
                        pass(ii)      = (sc == 1);
                        stabv(ii)     = (st == 1);
                        irf_list{ii}  = irf_tmp;     % store per-candidate
                        Qdraw_list{ii}  = Qdraw_tmp;     % store per-candidate
                        xrot_list{ii} = xrot_tmp;    % store per-candidate
                        lambda_list{ii} = lambda_tmp;    % store per-candidate
                        By_list{ii} = By_tmp;
                        fsign_list{ii} = fsign_tmp;
                    end
                    
                    % --- choose the first successful candidate
                    idx = find(pass & stabv, 1, 'first');
                    if ~isempty(idx)
                        % accept proposal
                        x0   = x1;
                        var0 = var1;
                        lh0  = lh1;                      % ok: MH accepted x1 earlier
                    
                        % if you also track these, make sure they exist BEFORE assigning here
                        likelihood0    = likelihood1;
                        A0_prior0      = A0_prior1;
                        lambda_prior0  = lambda_prior1;
                    
                        recorded_irfs = irf_list{idx};
                        Qdraw         = Qdraw_list{idx};
                        x1_rotated    = xrot_list{idx};
                        lambda        = lambda_list{idx};
                        By            = By_list{idx};
                        fsign         = fsign_list{idx};
                    
                        % Optional: recompute posterior at the rotated draw for diagnostics
                        [lh_wQ, likelihood_wQ, A0_prior_wQ, lambda_prior_wQ, ~] = eval_posterior( ...
                            x1_rotated, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
                    
                        acceptance_burning = acceptance_burning + 1;
                        SignCheck = 1;
                    else
                        SignCheck = 0;
                    end

                end
            end
    
            % === Gibbs step for delta ===
            eout = var0.udraw(1:size(delta0,1), :) .* exp(delta0);
            if isempty(tparam)
                delta0 = drawdelta(eout, alpha, K);
            else
                delta0 = drawt(eout, tparam, (tscale^2) * tparam);
            end
    
        % -------------------------------
        % Sign-restriction only, no heteroskedasticity
        % -------------------------------
        elseif use_signrestrictions && ~any(lcLmd(:))
    
            SignCheck = 0;
            nDrawAttempts = 0;
    
            while SignCheck == 0
                nDrawAttempts = nDrawAttempts + 1;
                x1 = x0 + mvnrnd(zeros(length(x0),1), Sigma)';
    
                [lh1, likelihood1, A0_prior1, lambda_prior1, var1] = eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);
                if isempty(var1) || ~isstruct(var1) || ~isfield(var1, 'Bdraw')
                    continue;
                end
    
                % --- pre-draw distinct Q's on the client (vary seed by iburn and/or proposal count)
                Ncand  = inneriteration;
                Qbatch = draw_Q_batch(ny, Ncand);   % your helper with optional seed
                
                % --- preallocate outputs that vary by ii
                pass      = false(Ncand,1);
                stabv     = false(Ncand,1);
                irf_list  = cell(Ncand,1);
                Qdraw_list = cell(Ncand,1);                
                xrot_list = cell(Ncand,1);
                lambda_list = cell(Ncand,1);
                By_list = cell(Ncand,1);
                fsign_list = cell(Ncand, 1);

                parfor ii = 1:Ncand
                    Qk = Qbatch(:,:,ii);   % deterministic, different per ii
                
                    [sc, ~, irf_tmp, ~, Qdraw_tmp, xrot_tmp, st, lambda_tmp, By_tmp, fsign_tmp] = SignRestrictionCheck( ...
                        x1, var1, ...
                        SignRestrictions, NarrativeRestrictions, ...
                        lags, lcA0, constant, ny, irf_horizon, ...
                        sign_regime_dependent, penalty_offdiagonal_on, lcLmd, fix_first_regime, Qk);
                
                    pass(ii)      = (sc == 1);
                    stabv(ii)     = (st == 1);
                    irf_list{ii}  = irf_tmp;     % store per-candidate
                    Qdraw_list{ii}  = Qdraw_tmp;     % store per-candidate
                    xrot_list{ii} = xrot_tmp;    % store per-candidate
                    lambda_list{ii} = lambda_tmp;    % store per-candidate
                    By_list{ii} = By_tmp;
                    fsign_list{ii} = fsign_tmp;
                end
                
                % --- choose the first successful candidate
                idx = find(pass & stabv, 1, 'first');
                if ~isempty(idx)
                    % accept proposal
                    x0   = x1;
                    var0 = var1;
                    lh0  = lh1;                      % ok: MH accepted x1 earlier
                
                    % if you also track these, make sure they exist BEFORE assigning here
                    likelihood0    = likelihood1;
                    A0_prior0      = A0_prior1;
                    lambda_prior0  = lambda_prior1;
                
                    recorded_irfs = irf_list{idx};
                    Qdraw         = Qdraw_list{idx};
                    x1_rotated    = xrot_list{idx};
                    lambda        = lambda_list{idx};
                    By            = By_list{idx};
                    fsign         = fsign_list{idx};

                    % Optional: recompute posterior at the rotated draw for diagnostics
                    [lh_wQ, likelihood_wQ, A0_prior_wQ, lambda_prior_wQ, ~] = eval_posterior( ...
                        x1_rotated, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
                
                    acceptance_burning = acceptance_burning + 1;
                    SignCheck = 1;
                else
                    SignCheck = 0;
                end               
            end

        % -------------------------------
        % No sign restrictions
        % -------------------------------
        else
            
            x1 = x0 + mvnrnd(zeros(length(x0), 1), Sigma)';
            if any(lcLmd(:))
                [lh1, likelihood1, A0_prior1, lambda_prior1, var1] = eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options);
            else
                [lh1, likelihood1, A0_prior1, lambda_prior1, var1] = eval_posterior(x1, y, lags, lcA0, lcLmd, Tsigbrk, [], tvA, options);
            end

            if log(rand()) < lh0 - lh1
                x0 = x1; var0 = var1; lh0 = lh1;
                likelihood0 = likelihood1; A0_prior0 = A0_prior1; lambda_prior0 = lambda_prior1;
            end

            if any(lcLmd(:))
                eout = var0.udraw(1:size(delta0, 1), :) .* exp(delta0);        
                if isempty(tparam)
                    delta0 = drawdelta(eout, alpha, K);
                else
                    delta0 = drawt(eout, tparam, (tscale^2) * tparam);
                end
            end
        end
    end
    acceptance_sampling = acceptance_sampling + 1;
    % fprintf("Accepted samplings are %d\n", acceptance_sampling);
    time_to_eachdraw = toc(samplein_start_time);
    fprintf('[INFO] Time to reach %4d accepted sampling-ins: %.2f seconds\n', idraw,time_to_eachdraw);

    % === Store ===
    if use_signrestrictions
        for k = 1:length(accepted_shock)
            shock_k = accepted_shock(k);
            sign_irfs_temp(:,:,shock_k,idraw) = recorded_irfs(:,:,shock_k);
        end
        draw_lh_wQ(idraw,:) = lh_wQ;
        draw_likelihood_wQ(idraw,:) = likelihood_wQ;
        draw_A0prior_wQ(idraw,:) = A0_prior_wQ;
        draw_lambdaprior_wQ(idraw,:) = lambda_prior_wQ;
        draw_Qdraw(:,:,idraw) = Qdraw;
        draw_Lambda(:,:,idraw) = lambda;
        draw_By(:,:,:,idraw) = By;
        draw_fsign(:,:,idraw) = fsign;
    end

    Phi = var0.Bdraw;
    draw_lh(idraw,:) = lh0;
    draw_likelihood(idraw,:) = likelihood0;
    draw_A0prior(idraw,:) = A0_prior0;
    draw_lambdaprior(idraw,:) = lambda_prior0;
    draw_x0(idraw,:) = x0;
    draw_Phi(:,:,idraw) = Phi;
    draw_udraw(:,:,idraw) = var0.udraw(1:T-lags,:);
    draw_udraw_true(:,:,idraw) = var0.udraw_true(1:T-lags,:);
    if any(lcLmd(:))
        draw_dout(:,:,idraw) = delta0;
    end
    
    if mod(idraw,10)==0
        fprintf('Sampling iteration: %4d / %4d | Acceptance: %.2f%%\n', ...
            idraw, ndraw, 100*acceptance_sampling/idraw);       
        time_to_10 = toc(samplein_start_time);
        fprintf('[INFO] Time to reach %4d accepted sampling-ins: %.2f seconds\n', idraw,time_to_10);
    end
end
% Final sampling summary
fprintf('Final acceptance number in sampling: %.5f\n', acceptance_sampling);

%% Draw IRFs by regime (scale-only vs. rotated-λ)
if nRegimes > 1
    xA_draws = draw_x0(:, 1:sum(lcA0(:)));               % ndraw x nAparams
    K = numel(accepted_shock);                           % constant length
    for i = 1:nRegimes
        irf_save = zeros(ny, irf_horizon, K, ndraw);  % preallocate for parfor

        parfor j = 1:ndraw
            % --- A0 draw ---
            A0 = zeros(ny);
            A0(lcA0) = xA_draws(j, :); 
            
            % ---By draw ---
            By = draw_By(:,:,:,j);

            % --- Q draw (columns are q_j) ---
            Q = draw_Qdraw(:, :, j);

            % --- regime variance ---
            Dhalf = diag(sqrt(max(draw_Lambda(:, i, j), 0)));  % Λ_r^{1/2}

            if sign_regime_dependent == 0
                % SCALE but do NOT rotate λ:  S = A0^{-1} Λ^{1/2} Q
                Smat = (A0 \ Q') * Dhalf;
            else
                % ROTATE λ explicitly: S = A0^{-1} Q (Q' Λ Q)^{1/2}
                % (This equals (A0 \ Dhalf) * Q, but computed explicitly if desired.)
                Smat = (A0 \ Q') * sqrtm(Q * diag(draw_Lambda(:, i, j)) * Q');
            end

            % --- IRFs for this draw ---
            irTemp = impulsdtrf(By, Smat, irf_horizon);

            % store only the accepted shocks, in compact [1..numel(accepted_shock)]
            for ij = 1:K
                irf_save(:, :, ij, j) = draw_fsign(1, ij, j) * irTemp(:, :, accepted_shock(ij));
            end
        end

        output.sign_irf_regime(:, :, :, :, i) = irf_save;
    end
end

%% Output
output.draw_lh = draw_lh;
output.draw_likelihood = draw_likelihood;
output.draw_A0prior = draw_A0prior;
output.draw_lambdaprior = draw_lambdaprior;
output.draw_x0 = draw_x0;
output.draw_Phi = draw_Phi;
output.draw_udraw = draw_udraw;
output.draw_udraw_true = draw_udraw_true;
output.lcA0 = lcA0;
output.lcLmd = lcLmd;
output.lrange = 1:nRegimes;
output.varnames = varnames;
if any(lcLmd(:))
    output.draw_dout = draw_dout;
end

if isfield(options, 'SignRestrictions') && ~isempty(options.SignRestrictions)
    output.sign_attempts = sign_attempts;
    output.sign_accepted_shock = sign_accepted_shock;
    output.sign_irfs_temp = sign_irfs_temp;
    output.sign_Qdraws = sign_Qdraws;
    output.draw_lh_wQ = draw_lh_wQ;
    output.draw_likelihood_wQ = draw_likelihood_wQ;
    output.draw_A0prior_wQ = draw_A0prior_wQ;
    output.draw_lambdaprior_wQ = draw_lambdaprior_wQ;
    output.draw_Qdraw = draw_Qdraw;
    output.draw_Lambda = draw_Lambda;
    output.draw_By = draw_By;
    output.draw_fsign = draw_fsign;
end

end

