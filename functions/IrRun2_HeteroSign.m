function output = IrRun2_HeteroSign(draw_x, y, lcA0, lcLmd, lags, nStep, lrange, Qdraw)
%--------------------------------------------------------------------------
% MATLAB translation and modification of IrRun2 function by Karthik Sastry
% This implement to calculate IRFs from Heteroskedasticity and Sign restri-
% ctions
% Inputs:
%   A, Aplus: structural matrices
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
startIdx = sum(lcA0(:)) + sum(lcLmd(:)) + 1;
AplusVec = draw_x(startIdx:startIdx + nvar*nvar*lags - 1);
Aplus = reshape(AplusVec, [nvar, nvar, lags]);

% Extract A0 matrix
A = zeros(nvar);
A(lcA0) = draw_x(1:sum(lcA0(:)));

% Extract lambda matrix
lambda = zeros(size(lcLmd));
lambda(lcLmd) = draw_x(sum(lcA0(:)) + (1:sum(lcLmd(:))));
lambda(:,end) = size(lambda,2) - sum(lambda(:,1:end-1),2);

ir = NaN(nvar, nStep, nvar, nRegimes);

% Loop over regimes (or dummy if no regime-switching)
for iL = 1:nRegimes
    % Direct indexing using sequential regime index
    Ainv = inv(A);

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
    if exist('Qdraw', 'var') && ~isempty(Qdraw) && all(~isnan(Qdraw(:)))
        smat_base = Ainv * Qdraw;
    else
        smat_base = Ainv;
    end
    
    irFull = impulsdtrf(tempvar, smat_base, nStep);
    ir(:,:,:,iL) = irFull;
end

output.ir = ir;
output.A0 = A(:,:,lrange);
output.Aplus = Aplus;
output.By = By;
output.lambda = lambda;
output.signShock = [];
end
