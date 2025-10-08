function output = IrRun2_tvA(draw_x, y, lcA0, lcLmd, lags, nStep, nA, lrange, SignRestrictions)
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

nvar = size(lcA0,1);
T = size(y,1);
nRegimes = length(lrange);

% Reshape A+ parameters
startIdx = nA + 1;
AplusVec = draw_x(startIdx:startIdx + nvar*nvar*lags - 1);
Aplus = reshape(AplusVec, [nvar, nvar, lags]);

% Extract A0 matrix
nA = nA -sum(lcLmd(:));
nregimes = length(draw_x(1:nA)) / sum(lcA0(:)); % Corrected
A = zeros(nvar, nvar, nregimes);
Acoef = reshape(draw_x(1:nA), sum(lcA0(:)), nregimes);

for r = 1:nregimes
    tempA = zeros(nvar, nvar);
    tempA(lcA0) = Acoef(:, r);
    A(:, :, r) = tempA;
end

% Extract lambda matrix
if size(lcLmd, 2) == 1 % No regimes case
    lambda = ones(nvar, 1);
else
    lambda = zeros(size(lcLmd));
    lambda(lcLmd) = draw_x(nA + (1:sum(lcLmd(:))));
    lambda(:,end) = size(lambda,2) - sum(lambda(:,1:end-1),2);
end

ir = NaN(nvar, nStep, nvar, nRegimes);
Ainv_all = NaN(nvar, nvar, nRegimes);

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
end

% Loop over regimes (or dummy if no regime-switching)
for iL = 1:nRegimes
    % Direct indexing using sequential regime index
    iLambda = lrange(iL);
    Ainv = inv(A(:,:,iLambda));
    Ainv_all(:,:,iL) = Ainv;

    % VAR coefficient matrices
    By = zeros(nvar, nvar, lags);
    for iLag = 1:lags
        By(:,:,iLag) = Ainv * Aplus(:,:,iLag);
    end
    tempvar.By = By;
    
    %Check if regime-switching volatility is present
    % if size(lcLmd, 2) == 1
    %     smat_base = Ainv;
    % else
    %     smat_base = Ainv* diag(sqrt(lambda(:, iLambda)));
    % end
    smat_base = Ainv;
    
    satisfied = false;
    maxTries = 10000;
    trial = 0;

    if ~useSign
        irFull = impulsdtrf(tempvar, smat_base, nStep);
        ir(:,:,:,iL) = irFull;
    else
        while ~satisfied && trial < maxTries
            trial = trial + 1;

            [Q, R] = qr(randn(nvar));
            for i = 1:nvar
                if R(i,i) < 0
                    Q(:,i) = -Q(:,i);
                end
            end

            smat_trial = smat_base * Q;
            irTemp = impulsdtrf(tempvar, smat_trial, signHorizon + lags);
            [d, accepted_shock, fsign] = checkrestrictions_allshocks(SignRestrictions, irTemp);
            if d == 1
                satisfied = true;
                irFull = impulsdtrf(tempvar, smat_trial, nStep);
                ir(:,:,:,iLambda) = fsign * irFull;
            end
        end

        if ~satisfied
            warning('Sign restrictions not satisfied after %d trials (regime %d).', maxTries, iLambda);
        end

        % % Try with smat_base as-is
        % irTemp = impulsdtrf(tempvar, smat_base, signHorizon + lags);
        % [d, fsign] = checkrestrictions2(signRestrictions, irTemp);
        % 
        % if d == 1
        %     irFull = impulsdtrf(tempvar, smat_base, nStep);
        %     ir(:,:,:,iL) = fsign * irFull;
        % else
        %     % Try with sign-flipped smat_base
        %     smat_trial = -1 * smat_base;
        %     irTemp = impulsdtrf(tempvar, smat_trial, signHorizon + lags);
        %     [d, fsign] = checkrestrictions2(signRestrictions, irTemp);
        % 
        %     if d == 1
        %         irFull = impulsdtrf(tempvar, smat_trial, nStep);
        %         ir(:,:,:,iL) = fsign * irFull;
        %     else
        %         warning('Sign restrictions not satisfied with smat_base or -smat_base (regime %d)', iLambda);
        %     end
        % end
    end
end

output.ir = ir;
output.A0 = A(:,:,lrange);
output.Ainv = Ainv_all;
output.Aplus = Aplus;
output.By = By;
%output.umat = umat;
output.lambda = lambda;

if useSign
    output.signShock = accepted_shock;
else
    output.signShock = [];
end

end
