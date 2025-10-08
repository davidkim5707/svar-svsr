function [log_posterior, log_likelihood, A0_prior_logprob, lplmd,var] = bvar_posterior(seedx, y, lags, lcA0, lcLmd, Tsigbrk,oweights, options)
% This MATLAB code is translated and modified from R code originally 
% written by Karthik Sastry, part of the SVAR toolkit.

% Settings
if nargin < 9
    oweights = [];
end

% SETTINGS
[T,ny] = size(y);
first_obs = lags + 1;
presample = 0;
lmdmean = 0;

%% Preparing A0 matrix
A = zeros(ny, ny);
nA = sum(lcA0(:));
A(lcA0) = seedx(1:nA);

% Prior on A0: Gaussian mean 0, var 100
A0_mean = 100;
A0_var = 200;  % because A0_var^2 = 100 → A0_var = 10

A0_prior_logprob = -0.5 * sum(((A - eye(ny) * A0_mean).^2) / (A0_var^2), 'all') ...
                   - ny^2 * (log(2 * pi)/2 + log(A0_var));

%% Preparing lmd
nLmd = sum(lcLmd(:));
nSig = length(Tsigbrk);

% Option flag (true/1 fixes first regime to 1)
fix_first = isfield(options,'fix_first_regime') && logical(options.fix_first_regime);

if nLmd > 0
    % Initialize
    lmd = ones(ny, nSig);  % start with ones so col-1=1 is already correct if fixed
    % Fill only free cells from seedx
    lmd(lcLmd) = seedx(nA + (1:nLmd));

    if fix_first
        % Hard-fix first regime to 1, never let params override it
        lmd(:,1) = 1;

        % Enforce sum over remaining regimes to be (nSig-1)
        if nSig >= 3
            partial = sum(lmd(:, 2:(nSig-1)), 2);                  % sum of middle free cols
            lmd(:, nSig) = (nSig - 1) - partial;                   % residual column
        else
            % nSig == 2: only one free regime left -> must be 1
            lmd(:, nSig) = 1;
        end

        % Positivity check (reject draw if violated)
        if any(lmd(:,2:end) <= 0, 'all')
            % => reject/penalize; avoid "max(...,eps)" because it biases the prior/likelihood
            disp('Warning: Non-positive lambda created with fix_first_regime; rejecting draw.');
            log_posterior = 1e5;          % your rejection convention
            log_likelihood = -1e10;
            lplmd = -1e5;
            var = [];
            return;
        end

        % Dirichlet(2) prior
        lpL = sum(log(lmd(:,2:end)) - log(nSig-1), 1) - gammaln(2 * (nSig-1));
        lplmd = sum(lpL) - (nSig - 1-1) * log(nSig-1);
    else
        % ===== ORIGINAL BEHAVIOR (no fixing of first regime) =====
        % Enforce arithmetic-average constraint so row-sum is nSig
        rowSums = sum(lmd(:, 1:(nSig-1)), 2);
        lmd(:, nSig) = nSig - rowSums;

        % Positivity check
        if any(lmd(:) <= 0)
            disp('Warning: Non-positive lambda; rejecting draw.');
            log_posterior = 1e5;
            log_likelihood = -1e10;
            lplmd = -1e5;
            var = [];
            return;
        end
        % Dirichlet(2) prior
        lpL = sum(log(lmd) - log(nSig), 1) - gammaln(2 * nSig);
        lplmd = sum(lpL) - (nSig - 1) * log(nSig);
    end
    

else
    lmd = ones(ny, nSig);
    lplmd = 0;
end

% Take log
llmd = -log(lmd);

% Append mean if needed
if isvector(llmd)
    lmdbar = llmd;
elseif lmdmean
    lmdbar = mean(llmd, 2);
else
    lmdbar = zeros(size(llmd,1),1);
end
llmd = [llmd, lmdbar];

sigpar.A0 = A;
sigpar.lmd = llmd;
sigpar.Tsigbrk = [Tsigbrk,T];

%% Hyperparameters
bvar_prior_tau = options.minn_prior_tau;
bvar_prior_decay = options.minn_prior_decay;
bvar_prior_lambda = options.minn_prior_lambda;
bvar_prior_mu = options.minn_prior_mu;
bvar_prior_omega = options.minn_prior_omega;

mnprior.tight = bvar_prior_tau;
mnprior.decay = bvar_prior_decay;
mnprior.unit_root_ = options.unitroot;

vprior.w = bvar_prior_omega;
if ~isempty(vprior.w)
    vprior.sig = std(y(first_obs-lags : first_obs+presample-1,:))';
else
    vprior.sig = [];
end

urprior.lambda = bvar_prior_lambda;
urprior.mu = bvar_prior_mu;

ybar = mean(y(1:lags, :), 1);

% Calculate the dimensions of y
[T, ~] = size(y);  % T is number of observations, nv is number of endogenous variables

% Add a constant column if options.constant == true
if options.constant
    xdata = ones(T, 1);  % T is the number of time observations
else
    xdata = [];
end

% Validate dimensions and compute nx (number of exogenous regressors)
if isempty(xdata)
    nx = 0;
else
    [Tx, nx] = size(xdata);
    if Tx ~= T
        error('Mismatch: xdata must have the same number of rows as ydata (T).');
    end
end

if isempty(mnprior.tight) && isempty(mnprior.decay) && isempty(mnprior.unit_root_) && ...
   isempty(urprior.lambda) && isempty(urprior.mu) && ...
   isempty(vprior.w)
    ydum = [];
    xdum = [];
    pbreaks = [];
else
    [ydum, xdum, pbreaks] = varprior(ny, nx, lags, mnprior, vprior,urprior,ybar);
end

%% Posterior likelihood calculation
ydata = y;

var = rfvar3([ydata; ydum], lags, [xdata; xdum], [T; T + pbreaks], [], [], sigpar, oweights);

% Log-likelihood
Tu = size(var.u,1);
lmdllh = 0.5 * sum(var.lmdseries, 'all');
detTerm = log(abs(det(A)));
llh = -0.5 * sum(var.u(:).^2) + Tu * (-ny * log(2 * pi) / 2 + detTerm) + lmdllh;
if options.constant
    nX = lags * ny + 1;
else
    nX = lags * ny;
end
log_dnsty = llh + 0.5 * sum(var.logdetxxi, 'all') + ny * nX * log(2 * pi) / 2;

%% Prior likelihood
Tp = presample + lags;
priorTsigbrk = [0, Tp];
priorlmd = [llmd(:,1), llmd(:, end)];

sigpar.lmd = priorlmd;
sigpar.Tsigbrk = priorTsigbrk;

if isempty(xdata)
    xprior = xdum;
else
    xprior = [xdata(1:Tp,:); xdum];
end

varp = rfvar3([ydata(1:Tp,:); ydum], lags, xprior, [Tp; Tp + pbreaks],  [], [], sigpar, []);

Tup = size(varp.u,1);
lmdllhp = 0.5 * sum(varp.lmdseries, 'all');
detPriorA0 = log(abs(det(A)));

llhp = -0.5 * sum(varp.u(:).^2) - Tup * (ny * log(2 * pi) / 2 - detPriorA0) + lmdllhp;
normalizer = 0.5 * sum(varp.logdetxxi, 'all') + ny * nX * log(2 * pi) / 2;
log_prior = llhp + normalizer;

%% Final posterior
log_likelihood = log_dnsty-log_prior;
log_posterior = -1*sum(log_likelihood, 'all') - A0_prior_logprob - lplmd;
end
