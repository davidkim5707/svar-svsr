function info = SetupZeroInfo(nvar, npredetermined, SignRestrictions)
% SETUPZEROINFO
%   Build zero-restriction bookkeeping in RRWZ style plus Z selection matrices.
%   W{j} has nvar - ((j-1) + r_j) rows, where r_j is the number of zero rows on shock j.
%
% Outputs (key fields)
%   info.zeros_per_shock : nvar x 1 (# zero horizons summed over variables, per shock j)
%   info.W{j}            : ( nvar - ((j-1) + zeros_per_shock(j)) ) x nvar
%   info.horizons        : vector (0-based) of horizons that appear in zero restrictions
%   info.details(k)      : .var, .shock, .tset (1-based t), .hset (0-based h), .raw
%   info.Z{j}            : (#zero rows for shock j) x ( numel(info.horizons)*nvar )
%   info.Z_var{j}        : vector (length = #zero rows for shock j), variable index for each row
%   info.Z_blockIdx{j}   : vector of block indices in info.horizons for each row

    if isstring(SignRestrictions), SignRestrictions = cellstr(SignRestrictions); end
    assert(iscell(SignRestrictions), 'SignRestrictions must be a cell array of strings.');

    % ---- Parse zero restrictions ----
    details = struct('var',{},'shock',{},'tset',[],'hset',[],'raw','');
    zeros_per_var   = zeros(nvar,1);
    zeros_per_shock = zeros(nvar,1);
    all_h = [];

    % regex: y(i, tset, s) = 0   where tset can be a, a:b, or a,b,c
    pat = 'y\(\s*(\d+)\s*,\s*([^,]+)\s*,\s*(\d+)\s*\)\s*=\s*0';

    for r = 1:numel(SignRestrictions)
        str = strtrim(SignRestrictions{r});
        if ~(contains(str,'=0') || contains(str,'= 0')), continue; end

        tok = regexp(str, pat, 'tokens', 'once');
        if isempty(tok), continue; end

        iVar   = str2double(tok{1});
        tSpec  = strtrim(tok{2});             % 1-based t in your DSL
        jShock = str2double(tok{3});

        tset = parseHorizonSpec(tSpec);       % e.g., [1] or 1:K or [1 3 5]
        if isempty(tset), continue; end
        if any(tset < 1)
            error('Horizon t must be >= 1 (got %s).', tSpec);
        end

        hset = tset - 1;                      % map DSL t -> IRF horizon h = t-1

        details(end+1) = struct( ...
            'var', iVar, 'shock', jShock, ...
            'tset', tset(:).', 'hset', hset(:).', 'raw', str); %#ok<AGROW>

        % counts (each horizon counts as one zero)
        if iVar>=1 && iVar<=nvar
            zeros_per_var(iVar) = zeros_per_var(iVar) + numel(hset);
        end
        if jShock>=1 && jShock<=nvar
            zeros_per_shock(jShock) = zeros_per_shock(jShock) + numel(hset);
        end
        all_h = [all_h, hset(:).']; %#ok<AGROW>
    end

    unique_h = unique(all_h);
    % If no zeros were specified, keep empty vector (caller can handle)
    % Otherwise, unique_h is 0-based horizons used by zeros.
    H = numel(unique_h);

    % ---- RRWZ sizing: s_j = nvar - ((j-1) + r_j) ----
    r = zeros_per_shock(:);                    % r_j
    s = nvar - ((0:nvar-1)' + r);              % s_j
    if any(s < 0)
        jbad = find(s < 0, 1);
        error('Infeasible zeros: column j=%d leaves free dim %d < 0 (r_j=%d).', jbad, s(jbad), r(jbad));
    end

    % ---- Build W (Gaussian rows) ----
    W = cell(nvar,1);
    for j = 1:nvar
        if s(j) == 0
            W{j} = zeros(0, nvar);
        else
            W{j} = randn(s(j), nvar);
        end
    end

    % ---- Build Z selection matrices (by shock) ----
    % Each Z{j} has one row per zero restriction (per horizon instance) on shock j.
    % Width = H * nvar, where H = numel(unique_h) (0-based horizons).
    Z = cell(nvar,1);
    Z_var = cell(nvar,1);        % variable index per row
    Z_blockIdx = cell(nvar,1);   % block index (which horizon in unique_h) per row
    for j = 1:nvar
        Z{j} = zeros(0, H*nvar);
        Z_var{j} = zeros(0,1);
        Z_blockIdx{j} = zeros(0,1);
    end

    % Helper: map horizon h to block index b in unique_h
    % (b runs 1..H, IRF rows for shock j at horizon h sit at row index (b-1)*nvar + j)
    for k = 1:numel(details)
        iVar   = details(k).var;
        jShock = details(k).shock;
        hset   = details(k).hset;

        for hh = hset
            if isempty(unique_h)
                error('No horizons accumulated for zeros, but a zero was parsed.'); % defensive
            end
            b = find(unique_h == hh, 1, 'first');
            if isempty(b)
                error('Internal: horizon h=%d not found in info.horizons.', hh);
            end

            zrow = zeros(1, H*nvar);
            row_in_stack = (b-1)*nvar + jShock;   % selects the stacked-IRF row
            zrow(1, row_in_stack) = 1;

            Z{jShock}        = [Z{jShock}; zrow];           
            Z_var{jShock}    = [Z_var{jShock}; iVar];       
            Z_blockIdx{jShock}= [Z_blockIdx{jShock}; b];    
        end
    end

    % ---- Pack info ----
    info.nvar              = nvar;
    info.npredetermined    = npredetermined;
    info.nzeros            = sum(r);
    info.dim               = sum(s);
    info.horizons          = unique_h;          % 0-based horizons used by zeros
    info.zeros_per_var     = zeros_per_var;
    info.zeros_per_shock   = r;
    info.details           = details;           % includes both 1-based t and 0-based h
    info.W                 = W;

    % Z-related outputs
    info.Z                 = Z;                 % selection matrices
    info.Z_var             = Z_var;             % variable index per row of Z{j}
    info.Z_blockIdx        = Z_blockIdx;        % which horizon-block each row targets
end

% --------- Helper: parse horizon specs like '1', '1:3', '1,3,5' ----------
function h = parseHorizonSpec(hSpec)
    hs = regexprep(hSpec, '\s+', '');
    if isempty(hs), h = []; return; end
    if contains(lower(hs),'inf')
        error('Zero restrictions with t = Inf are not supported.');
    end
    if contains(hs, ':')
        c = str2double(split(hs,':'));
        if numel(c)==2 && all(isfinite(c))
            h = c(1):c(2);
        else
            h = []; return;
        end
    elseif contains(hs, ',')
        h = str2double(split(hs, ','));
    else
        h = str2double(hs);
    end
    h = unique(round(h(:).'));
    if any(isnan(h)), h = []; end
end
