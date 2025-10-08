function var  = rfvar3(ydata, lags, xdata, breaks, lambda, mu, sigpar, oweights)
%function var=rfvar3(ydata,lags,xdata,breaks,lambda,mu)
% This algorithm goes for accuracy without worrying about memory requirements.
% ydata:   dependent variable data matrix
% xdata:   exogenous variable data matrix
% lags:    number of lags
% breaks:  rows in ydata and xdata after which there is a break.  This allows for
%          discontinuities in the data (e.g. war years) and for the possibility of
%          adding dummy observations to implement a prior.  This must be a column vector.
%          Note that a single dummy observation becomes lags+1 rows of the data matrix,
%          with a break separating it from the rest of the data.  The function treats the 
%          first lags observations at the top and after each "break" in ydata and xdata as
%          initial conditions. 
% lambda:  weight on "co-persistence" prior dummy observations.  This expresses
%          belief that when data on *all* y's are stable at their initial levels, they will
%          tend to persist at that level.  lambda=5 is a reasonable first try.  With lambda<0,
%          constant term is not included in the dummy observation, so that stationary models
%          with means equal to initial ybar do not fit the prior mean.  With lambda>0, the prior
%          implies that large constants are unlikely if unit roots are present.
% mu:      weight on "own persistence" prior dummy observation.  Expresses belief
%          that when y_i has been stable at its initial level, it will tend to persist
%          at that level, regardless of the values of other variables.  There is
%          one of these for each variable.  A reasonable first guess is mu=2.
%      The program assumes that the first lags rows of ydata and xdata are real data, not dummies.
%      Dummy observations should go at the end, if any.  If pre-sample x's are not available,
%      repeating the initial xdata(lags+1,:) row or copying xdata(lags+1:2*lags,:) into 
%      xdata(1:lags,:) are reasonable subsititutes.  These values are used in forming the
%      persistence priors.
% A0:       (Optional) Structural identification matrix (Leeper & Roush, 2003).

% Original file downloaded from:
% http://sims.princeton.edu/yftp/VARtools/matlab/rfvar3.m

% Settings
if nargin < 8 
    oweights = [];
end

[T,nvar] = size(ydata);

nox = isempty(xdata);
if ~nox
    [T2,nx] = size(xdata);
else
    T2 = T;
    nx = 0;
    xdata = zeros(T2,0);
end

% note that x must be same length as y, even though first part of x will not be used.
% This is so that the lags parameter can be changed without reshaping the xdata matrix.
if T2 ~= T, error('Mismatch of x and y data lengths'),end
if nargin < 4
    nbreaks = 0;
    breaks = [];
else
    nbreaks = length(breaks);
end
breaks = [0;breaks;T];
smpl = [];
for nb = 1:nbreaks+1
    smpl = [smpl;[breaks(nb)+lags+1:breaks(nb+1)]'];
end
Tsmpl = size(smpl,1);
X = zeros(Tsmpl,nvar,lags);
for is = 1:length(smpl)
    X(is,:,:) = ydata(smpl(is)-(1:lags),:)';
end
X = [X(:,:) xdata(smpl,:)];
y = ydata(smpl,:);

% Everything now set up with input data for y=Xb+e 
% Add persistence dummies
if (~isempty(lambda) && lambda ~= 0) || (~isempty(mu) && mu > 0)
    ybar = mean(ydata(1:lags,:),1);
    if ~nox
        xbar = mean(xdata(1:lags,:),1);
    else
        xbar = [];
    end
    if lambda ~= 0
        if lambda>0
            xdum = lambda*[repmat(ybar,1,lags) xbar];
        else
            lambda = -lambda;
            xdum = lambda*[repmat(ybar,1,lags) zeros(size(xbar))];
        end
        ydum = zeros(1,nvar);
        ydum(1,:) = lambda*ybar;
        y = [y;ydum];
        X = [X;xdum];
    end
    if mu>0
        xdum = [repmat(diag(ybar),1,lags) zeros(nvar,nx)]*mu;
        ydum = mu*diag(ybar);
        X = [X;xdum];
        y = [y;ydum];
    end
end

if exist('sigpar', 'var') && ~isempty(sigpar)
    A0 = sigpar.A0;
    lmd = sigpar.lmd;
    Tsigbrk = sigpar.Tsigbrk;

    % Assume Tsigbrk is either empty ([]) or a vector of break dates (time indices)
    if ~isempty(Tsigbrk)
        % Tsigbrk is already time indices; no need for invtime()
        nsig = length(Tsigbrk);
    else
        nsig = 1;
    end
    
    % Finalize regime breakpoints
    Tsigbrk = [Tsigbrk, T];
    
    % Build regime index vector lmdndx
    lmdndx = [];
    for i = 1:nsig
        lmdndx = [lmdndx, repmat(i, 1, Tsigbrk(i+1) - Tsigbrk(i))];
    end
    
    % Now pad lmdndx if needed
    if max(smpl) > length(lmdndx)
        lastregime = lmdndx(end);
        lmdndx = [lmdndx, repmat(lastregime, 1, max(smpl) - length(lmdndx))];
    end
    
    % Safe assignment
    lmdseries = lmd(:, lmdndx);
    
    % Check if dummy observations were used in y (i.e., dim(y,1) > Tsmpl)
    if Tsmpl < size(y, 1)
        % Dummy observations included, so need to expand lmdseries to match full y size
        lmdp = mean(lmdseries(:, smpl), 2);  % mean across sample time points (smpl)
        % Extend lmdseries with mean values for dummy periods
        lmdseries = [lmdseries(:, smpl), repmat(lmdp, 1, size(y,1) - Tsmpl)];
    else
        % Just restrict lmdseries to sample points
        lmdseries = lmdseries(:, smpl);
    end

    if ~isempty(oweights)
        % Display info message
        %fprintf('oweights are applied in rfvar3\n');
        
        % Show oweights
        %fprintf('oweights:\n');
        %disp(oweights(1:3,1:3));
        
        % Show lmdseries BEFORE applying oweights
        colrange = 1:size(oweights, 1);
        %fprintf('lmdseries before applying oweights (affected columns):\n');
        %disp(lmdseries(1:3, 1:3));
        
        % Apply the transformation (log std dev scaling)
        lmdseries(:, colrange) = lmdseries(:, colrange) - 2 * oweights';
        
        % Show lmdseries AFTER applying oweights
        %fprintf('lmdseries after applying oweights (affected columns):\n');
        %disp(lmdseries(1:3, 1:3));
    end

    if ismatrix(A0)
        ya0 = y * A0';  % constant A0: same transformation for all
    else
        % Time-varying A0 case
        A0list = A0(:,:,lmdndx);   % select regime-specific A0s for each time
        A0list = A0list(:,:,smpl); % restrict to sample time points
    
        ya0 = zeros(size(y));
        for iy = 1:size(y,1)
            ya0(iy,:) = y(iy,:) * A0list(:,:,iy)';
        end
    end

    freqs = lmdndx(smpl);  % regime index for each observation in the sample

    listOutput = cell(1, nvar);
    errorflags = false(1, nvar);
    
    % Serial
    %tic
    for iq = 1:nvar
        [errorflag, coefs, u, logdetxxi, snglty, xx, coefs_draw, udraw, udraw_true] = ...
            linreg(iq, lmdseries, X, ya0, 1);
        errorflags(iq) = errorflag;
        listOutput{iq} = struct('coefs',coefs,'u',u,'logdetxxi',logdetxxi,'snglty',snglty,'xx',xx, ...
                                'coefs_draw',coefs_draw,'udraw',udraw,'udraw_true',udraw_true);
    end
    %t_serial = toc;
    
    % % Parallel (only for comparison)
    % pp = gcp('nocreate'); if isempty(pp), parpool; end
    % tic
    % parfor iq = 1:nvar
    %     [errorflag, coefs, u, logdetxxi, snglty, xx, coefs_draw, udraw, udraw_true] = ...
    %         linreg(iq, lmdseries, X, ya0, 1);
    %     errorflags(iq) = errorflag;
    %     listOutput{iq} = struct('coefs',coefs,'u',u,'logdetxxi',logdetxxi,'snglty',snglty,'xx',xx, ...
    %                             'coefs_draw',coefs_draw,'udraw',udraw,'udraw_true',udraw_true);
    % end
    % t_par = toc;
    % 
    % fprintf('serial=%.3fs  |  parfor=%.3fs\n', t_serial, t_par);
    
    if any(errorflags)
        warning('Error in linreg at equation(s): %s', mat2str(find(errorflags)));
        var = [];  % Empty output to prevent MATLAB error
        return;
    end
    
    % Preallocate
    ncoef = size(listOutput{1}.coefs, 1);
    nobs_u = size(listOutput{1}.u, 1);
    nX = size(listOutput{1}.xx, 1);
    
    B      = zeros(ncoef, nvar);
    u      = zeros(nobs_u, nvar);
    logdetxxi = NaN(nvar, 1);  % Use NaN to highlight skipped draws
    snglty = zeros(nvar, 1);
    Bdraw  = NaN(ncoef, nvar);
    udraw  = NaN(nobs_u, nvar);
    udraw_true  = NaN(nobs_u, nvar);
    xx     = NaN(nX, nX, nvar);
    
    for iv = 1:nvar
        B(:, iv) = listOutput{iv}.coefs;
        u(:, iv) = listOutput{iv}.u;
        xx(:, :, iv) = listOutput{iv}.xx;
    
        % Handle optional fields safely
        if isfield(listOutput{iv}, 'logdetxxi') && ~isempty(listOutput{iv}.logdetxxi)
            logdetxxi(iv) = listOutput{iv}.logdetxxi;
        end
    
        if isfield(listOutput{iv}, 'snglty')
            snglty(iv) = listOutput{iv}.snglty;
        end
    
        if isfield(listOutput{iv}, 'coefs_draw') && ~isempty(listOutput{iv}.coefs_draw)
            Bdraw(:, iv) = listOutput{iv}.coefs_draw;
        end
        
        if isfield(listOutput{iv}, 'udraw') && ~isempty(listOutput{iv}.udraw)
            udraw(:, iv) = listOutput{iv}.udraw;
        end

        if isfield(listOutput{iv}, 'udraw_true') && ~isempty(listOutput{iv}.udraw_true)
            udraw_true(:, iv) = listOutput{iv}.udraw_true;
        end
    end
    
else

    % Compute OLS regression and residuals
    [vl,d,vr] = svd(X,0);
    di = 1./diag(d);
    B = (vr.*repmat(di',nvar*lags+nx,1))*vl'*y;
    u = y-X*B;
    xxi = vr.*repmat(di',nvar*lags+nx,1);
    xxi = xxi*xxi';
    logdetxxi = 2 * sum(log(abs(di)));  % Compute log determinant
    
end

% Assign common fields
var.B = B;
var.u = u;
var.logdetxxi = logdetxxi;

if ~isempty(sigpar)  % Equivalent to !is.null(sigpar)
    var.udraw = udraw;
    var.udraw_true = udraw_true;
    var.xx = xx;
    var.lmdseries = lmdseries;
    var.Bdraw = Bdraw;
    var.udraw = udraw;
    var.snglty = snglty;
    var.freqs = freqs;
else
    var.xxi = xxi;
end
