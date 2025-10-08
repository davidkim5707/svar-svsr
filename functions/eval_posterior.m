function [lh, likelihood, A0_prior, lambda_prior, var] = eval_posterior(x, y, lags, lcA0, lcLmd, Tsigbrk, delta0, tvA, options)
%==========================================================================
% Evaluate the log posterior and associated output structure for a BVAR model.
%
% Inputs:
%   - x        : Current draw of model parameters (vectorized)
%   - y        : Observed data matrix (T x n)
%   - lags     : Number of VAR lags
%   - lcA0     : Logical matrix for identifying free elements in A0
%   - lcLmd    : Logical matrix for identifying regime-specific variances
%   - Tsigbrk  : Vector of regime breakpoints for stochastic volatility
%   - delta0   : Current draw of log-volatility parameters (used if lcLmd is nonzero)
%   - tvA      : Boolean flag, true if A0 matrix is time-varying
%
% Outputs:
%   - lh       : Log-posterior (up to a constant) of the current draw
%   - var      : Structure containing IRFs, residuals, fitted values, etc.
%
% This function selects the appropriate posterior evaluation routine depending
% on the heteroskedasticity structure (lcLmd) and whether A0 is time-varying (tvA).
%==========================================================================
if tvA
    if any(lcLmd(:))
        [lh, likelihood, A0_prior, lambda_prior, var] = bvar_posterior_tvA(x, y, lags, lcA0, lcLmd, Tsigbrk, delta0);
    else
        [lh, likelihood, A0_prior, lambda_prior, var] = bvar_posterior_tvA(x, y, lags, lcA0, lcLmd, Tsigbrk, options);
    end
elseif any(lcLmd(:))
    [lh, likelihood, A0_prior, lambda_prior, var] = bvar_posterior(x, y, lags, lcA0, lcLmd, Tsigbrk, delta0, options);
else
    [lh, likelihood, A0_prior, lambda_prior, var] = bvar_posterior(x, y, lags, lcA0, lcLmd, Tsigbrk, [], options);
end

end
