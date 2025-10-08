function output = IrRun2(draw_x, y, lcA0, lcLmd, lags, nStep, lrange, SignRestrictions)
%--------------------------------------------------------------------------
% MATLAB translation and modification of IrRun2 function by Karthik Sastry
% with sign restrictions integration inspired by bvar_ toolkit approach.
% Includes random rotation matrix (Q) following Rubio-Ramirez (2010).
% Also handles no-regime case using xx matrix for Cholesky IRFs.
% Inputs:
%   A, Aplus: structural matrices
%   umat: reduced-form residuals
%   lambda: diagonal elements for variance scaling
%   lrange: regime indices
%   lcLmd: logical matrix for regime-switching lambda
%   lags: number of lags
%   nStep: IRF horizon
%   signRestrictions: cell array of sign restriction structures (optional)
%--------------------------------------------------------------------------

nVar = size(lcA0,1);
T = size(y,1);
nRegimes = length(lrange);

% Reshape A+ parameters
startIdx = sum(lcA0(:)) + sum(lcLmd(:)) + 1;
AplusVec = draw_x(startIdx:startIdx + nVar*nVar*lags - 1);
Aplus = reshape(AplusVec, [nVar, nVar, lags]);

% Extract A0 matrix
A = zeros(nVar);
A(lcA0) = draw_x(1:sum(lcA0(:)));

% Extract lambda matrix
if size(lcLmd, 2) == 1 % No regimes case
    lambda = ones(nVar, 1);
else
    lambda = zeros(size(lcLmd));
    lambda(lcLmd) = draw_x(sum(lcA0(:)) + (1:sum(lcLmd(:))));
    lambda(:,end) = size(lambda,2) - sum(lambda(:,1:end-1),2);
end

ir = NaN(nVar, nStep, nVar, length(lrange));
Ainv_all = NaN(nVar, nVar, length(lrange));

if ~isempty(SignRestrictions)
    useSign = true;

    signHorizon = 0;
    for ii = 1:length(SignRestrictions)
        restr = SignRestrictions{ii};

        % Extract horizon from pattern like 1:h
        tokens = regexp(restr, '1:(\d+)', 'tokens');
        if ~isempty(tokens)
            horizonEnd = str2double(tokens{1}{1});
            signHorizon = max(signHorizon, horizonEnd);
        else
            signHorizon = max(signHorizon, nStep);
        end
    end
else
    useSign = false;
    signShock = NaN;  % or 0, if you prefer
end

% Precompute inverse of A
Ainv = inv(A);

% Loop over regimes (or dummy if no regime-switching)
for iLambda = 1:nRegimes
    regime = lrange(iLambda);
    Ainv_all(:,:,regime) = Ainv;

    % Prepare VAR coefficient matrices
    By = zeros(nVar, nVar, lags);
    for iLag = 1:lags
        By(:,:,iLag) = Ainv * Aplus(:,:,iLag);
    end
    tempvar.By = By;

    smat_base = Ainv;  % No heteroskedasticity
    tempvar.u = smat_base;

    if ~useSign
        % No sign restrictions → standard Cholesky IRFs
        ir_chol = impulsdtrf(tempvar, smat_base, nStep);
        ir(:,:,:,iLambda) = ir_chol;
    else

        % Sign restrictions active
        satisfied = false;
        maxTries = 10000;
        trial = 0;

        while ~satisfied && trial < maxTries
            trial = trial + 1;

            % Draw random orthonormal Q matrix
            [Q, R] = qr(randn(nVar));
            for i = 1:nVar
                if R(i,i) < 0
                    Q(:,i) = -Q(:,i);
                end
            end

            smat_trial = smat_base * Q;

            % Compute temporary IRF for checking sign restrictions
            extraSteps = 1;
            irTemp = impulsdtrf(tempvar, smat_trial, signHorizon + lags - 1 + extraSteps);
            [d, accepted_shock, fsign] = checkrestrictions_allshocks(SignRestrictions, irTemp);
            if d == 1
                satisfied = true;
                irFull = impulsdtrf(tempvar, smat_trial, nStep);
                ir(:,:,:,iLambda) = fsign * irFull;
            end
        end

        if ~satisfied
            warning('Sign restrictions not satisfied after %d trials (no-switching case).', maxTries);
        end
    end
end


output.ir = ir;
output.A0 = A;
output.Aplus = Aplus;
%output.umat = umat;
output.lambda = lambda;
output.Ainv = Ainv_all;
output.By = By;

if useSign
    output.signShock = accepted_shock;
else
    output.signShock = [];
end
end
