function [errorflag, coefs, u, logdetxxi, snglty, xx, coefs_draw, udraw, udraw_true] = ...
    linreg(iq, lmdseries, X, ya0, drawbe)
% This MATLAB code is translated and modified from R code originally 
% written by Karthik Sastry, part of the SVAR toolkit.

errorflag = false;

wt = exp(0.5 * lmdseries(iq, :));
wt = max(wt, 1e-10); % Increased safeguard for small weights
Xq = bsxfun(@times, X, wt');  % Weighted X
yq = wt' .* ya0(:, iq);       % Weighted y

% Check for bad weights and data
if any(isinf(Xq(:)))
    disp('Xq contains Inf values!');
    disp(find(isinf(Xq)));  % Print the indices
end

if any(isnan(Xq(:)))
    disp('Xq contains NaN values!');
    disp(find(isnan(Xq)));  % Print the indices
end

if any(isinf(yq))
    disp('yq contains Inf values!');
    disp(find(isinf(yq)));  % Print indices
end

if any(isnan(yq))
    disp('yq contains NaN values!');
    disp(find(isnan(yq)));  % Print indices
end

if any(wt < eps)
    disp('Weights contain values smaller than eps!');
    disp(find(wt < eps));  % Print indices
end

% If ANY bad condition triggers, set errorflag
if any(isinf(Xq(:))) || any(isnan(Xq(:))) || any(isinf(yq)) || any(isnan(yq)) || any(wt < eps)
    disp('Setting errorflag TRUE due to bad weights/data.');
    errorflag = true;
    coefs = []; u = []; logdetxxi = []; snglty = []; xx = []; coefs_draw = []; udraw = [];
    return;
end

% --- OLS estimation with Bayesian draw support and stability checks ---

% OLS
[Q, R] = qr(Xq, 0);
coefs = (R + 1e-6 * eye(size(R))) \ (Q' * yq); % Increased ridge regularization
u = yq - Xq * coefs;
xx = R' * R + 1e-6 * eye(size(R,2)); % Increased ridge regularization

% Draw coefficients
if drawbe
    [V, D] = eig(xx);
    eigvals = diag(D);
    if any(eigvals <= 1e-14)
        % Ill-conditioned → skip draw
        coefs_draw = []; udraw = []; snglty = true; logdetxxi = -Inf;
    else
        % Sample using inverse of sqrt(xx)
        sqrt_inv_xx = V * diag(1 ./ sqrt(eigvals)) * V';
        coefs_noise = sqrt_inv_xx * randn(size(xx, 1), 1);
        coefs_draw = coefs + coefs_noise;
        udraw = u - Xq * coefs_noise;
        udraw_true = udraw./wt';
        logdetxxi = -sum(log(eigvals));
        snglty = false;
    end
else
    coefs_draw = NaN(size(coefs));
    udraw = NaN(size(u));
    udraw_true = udraw/wt';
    logdetxxi = -2 * sum(log(abs(diag(R))));
    snglty = (logdetxxi == -Inf);
end

end
