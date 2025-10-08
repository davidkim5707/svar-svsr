function [log_posterior, var] = bvar_posterior_tvA(seedx, y, lags, lcA0, lcLmd, Tsigbrk,oweights)
% This MATLAB code is translated and modified from R code originally 
% written by Karthik Sastry, part of the SVAR toolkit.

% Settings
if nargin < 9
    oweights = [];
end

% SETTINGS
[T, ny] = size(y);
first_obs = lags + 1;
presample = 0;
lmdmean = 0;

%% Number of regimes
nSig = length(Tsigbrk);  % analogous to R's nSig
nA = sum(lcA0(:));       % number of non-fixed parameters per regime

%% Preparing A0 matrices across regimes
A = zeros(ny, ny, nSig);
Acoef = reshape(seedx(1:(nA * nSig)), nA, nSig);  % extract regime-specific A0 parameters

for iSig = 1:nSig
    Amat = zeros(ny, ny);
    Amat(lcA0) = Acoef(:, iSig);  % fill free parameters
    A(:, :, iSig) = Amat;
end

%% Prior on A0: Gaussian prior, mean 100, variance 200
allh = zeros(nSig, 1);
for iSig = 1:nSig
    Ascore = ((A(:, :, iSig) - 100 * eye(ny)).^2) / 4e4;
    allh(iSig) = -0.5 * sum(Ascore(lcA0)) - 0.5 * nA * (log(2 * pi) + log(200));
end
A0_prior_logprob = sum(allh);  % total log prior across regimes

%% Preparing lmd
nLmd = sum(lcLmd(:));
nSig = length(Tsigbrk);
nx = length(seedx);
hasLmd = (nx > nA * nSig);

if hasLmd
    % Fill in lambda parameters from the rest of seedx
    lmd = zeros(ny, nSig);
    lmd(lcLmd) = seedx(nA * nSig + (1:nLmd));
    
    % Enforce last regime's lmd column so rows sum to nSig
    rowSums = sum(lmd(:, 1:(nSig - 1)), 2);
    lmd(:, nSig) = nSig - rowSums;

    % Stability check: no negative variances
    if any(lmd(:) < 0)
        disp('Warning: Negative variance detected in lmd. Skipping this draw.');
        disp(['Min lmd value: ', num2str(min(lmd(:)))]);
        log_posterior = 1e5;
        var = []; % ensure var is assigned!
        return;
    end

    % Dirichlet(2) prior for each row
    % Prior: log(p(lmd)) ∝ sum(log(lmd)) + const
    lpL = sum(log(lmd) - log(nSig), 2) - gammaln(2 * nSig);  % row-wise
    lplmd = sum(lpL) - (nSig - 1) * log(nSig);               % total
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

if ndims(A) == 3
    A_mean = mean(A, 3);                     % average across regimes
    A = cat(3, A, A_mean);                  % append mean slice as new regime
end

sigpar.A0 = A;
sigpar.lmd = llmd;
sigpar.Tsigbrk = [Tsigbrk,T];

%% Hyperparameters
bvar_prior_tau = 3;
bvar_prior_decay = 0.5;
bvar_prior_lambda = 5;
bvar_prior_mu = 1;
bvar_prior_omega = 0;

mnprior.tight = bvar_prior_tau;
mnprior.decay = bvar_prior_decay;
mnprior.unit_root_ = ones(ny,1);

vprior.sig = std(y(first_obs-lags : first_obs+presample-1,:))';
vprior.w = bvar_prior_omega;

urprior.lambda = bvar_prior_lambda;
urprior.mu = bvar_prior_mu;

ybar = mean(y(1:lags, :), 1);

% Calculate the dimensions of y
[T, ~] = size(y);  % T is number of observations, nv is number of endogenous variables

% Add a constant term if const = true
xdata = ones(T, 1);

% Check if xdata is empty or not and calculate nx
if isempty(xdata)
    nx = 0;  % No exogenous variables if xdata is empty
else
    [Tx, nx] = size(xdata);  % Tx is the number of observations in xdata, nx is number of exogenous variables
end

% Safety check to ensure xdata length matches y
if ~isempty(xdata) && Tx ~= T
    error('The number of observations in xdata and y must match.');
end

[ydum, xdum, pbreaks] = varprior(ny, nx, lags, mnprior, vprior,urprior,ybar);

%% Posterior likelihood calculation
ydata = y;

var = rfvar3([ydata; ydum], lags, [xdata; xdum], [T; T + pbreaks], [], [], sigpar, oweights);

% Log-likelihood
Tu = size(var.u,1);
lmdllh = 0.5 * sum(var.lmdseries(:));

if ismatrix(A)
    detTerm = log(abs(det(A)));
else
    nSig = size(A, 3);
    dets = zeros(nSig, 1);
    for iSig = 1:nSig
        dets(iSig) = log(abs(det(A(:,:,iSig))));
    end

    % Use regime index per obs to select the correct determinant
    detTerm = mean(dets(var.freqs));  % var.freqs = regime index for each obs
end

llh = -0.5 * sum(var.u(:).^2) + Tu * (-ny * log(2 * pi) / 2 + detTerm) + lmdllh;
nX = lags * ny + 1;
log_dnsty = llh + 0.5 * sum(var.logdetxxi) + ny * nX * log(2 * pi) / 2;

%% Prior likelihood
Tp = presample + lags;
% Time-varying variances: set dummy breakpoints for prior structure
priorTsigbrk = [0, Tp];

% Prior A0 matrix: take first and last regime if regime-varying
if ndims(A) == 3
    priorA0 = cat(3, A(:,:,1), A(:,:,end));
else
    priorA0 = A;
end

% Prior lambda (variance scale): first and last regime's log lambda
priorlmd = [llmd(:,1), llmd(:,end)];

% Package prior structure
sigpar.A0 = priorA0;
sigpar.lmd = priorlmd;
sigpar.Tsigbrk = priorTsigbrk;

varp = rfvar3([ydata(1:Tp,:); ydum], lags, [xdata(1:Tp,1); xdum], [Tp; Tp + pbreaks],  [], [], sigpar, []);
Tup = size(varp.u,1);
lmdllhp = 0.5 * sum(varp.lmdseries(:));
if ndims(A) == 3
    % Only first and last regime used in priorA0
    detsp = [log(abs(det(A(:,:,1)))); log(abs(det(A(:,:,end))))];

    % Each obs in freqs ∈ {1,2}, corresponding to A0(:,:,1) or A0(:,:,end)
    % So we map regime index to detsp via freqs
    detPriorA0 = mean(detsp(varp.freqs));  % like detsp[varp$freqs] in R
else
    detPriorA0 = log(abs(det(A)));
end

llhp = -0.5 * sum(varp.u(:).^2) - Tup * (ny * log(2 * pi) / 2 - detPriorA0) + lmdllhp;
normalizer = 0.5 * sum(varp.logdetxxi) + ny * (lags * ny + 1) * log(2 * pi) / 2;
log_prior = llhp + normalizer;

%% Final posterior
log_posterior = log_dnsty-log_prior;
log_posterior = -1*log_posterior - A0_prior_logprob - lplmd;

end
