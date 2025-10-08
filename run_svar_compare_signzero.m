%********************************************************
%% Data loading and setting
%********************************************************
clc;        % Clear command window
clear;      % Clear all variables from workspace
close all;  % Close all figures

% Add the "functions" folder to MATLAB path
addpath(genpath('functions'));

% Load data from CSV file
filename = 'dataset_RubioRamirez_SignZero';
data = readtable(['data/' filename '.csv']);

% Extract variable names from column 2 to the last column
varnames = data.Properties.VariableNames(3:end);

% Convert date column to datetime format
decYears = data.date;        % e.g. 1955.1, 1955.2, ...
years  = floor(decYears);    % whole year part
qtr    = round((decYears - years)*10);  % 1,2,3,4

% Map quarter to first month
months = (qtr-1)*3 + 1;

% Construct datetime
data.dates = datetime(years, months, 1);

% Filter data to keep only observations within the date range
filtered_data = data(:, [end, 3:end-1]);
filtered_data.dates = datetime(filtered_data.dates, 'InputFormat', 'yyyy-MM-dd');

% Drop the 'date' column and define y
y = table2array(filtered_data(:,2:end));

%********************************************************
%% Define Options Struct for SVAR – Kilian-style Model
%********************************************************
% === USER SETTINGS ===
% --------------------------------------------------------------
% Uhlig (2005):             regime_on = 0, prior_on = 0, svar_svsr_on = 0
% Carriero et al. (2024):   regime_on = 1, prior_on = 0, svar_svsr_on = 0
% Dawis' method:            regime_on = 1, prior_on = 0, svar_svsr_on = 1
% --------------------------------------------------------------
regime_on     = 0;    % 1 = enable regime-switching
prior_on      = 0;    % 1 = use Minnesota prior
svar_svsr_on  = 0;    % 1 = Dawis (rotational invariance)
penalty_offdiagonal_on = 0; % 1 = Use penalty function to check sign-restrictions

% === THREADS SETTINGS ===
% if isempty(gcp('nocreate'))
%     parpool(8);  % Only create if none exists
% elseif gcp().NumWorkers ~= 8
%     delete(gcp('nocreate'));
%     parpool(8);
% end

% === BASIC MODEL SETTINGS ===
options.ny            = size(y, 2);    % Number of observable variables
options.lags          = 4;            % Number of lags
options.irf_horizon   = 40;            % IRF horizon
options.constant      = 1;             % Include constant
options.timetrend     = 0;             % No time trend
options.firstobs      = options.lags + 1;

% === SAMPLING AND COMPUTATION SETTINGS ===
options.presample     = 0;
options.ndraw         = 1000;         % Number of MCMC draws
options.nsep          = 1;            % Save every draw
options.nburn         = 1000;         % Burn-in draws
options.max_compute   = 1;            % Optimization: 1 = csminwel
options.nit           = 50;           % Iterations
options.hsnscale      = 0.01;          % Hessian scaling

% === DATA DIMENSIONS ===
options.nobs          = size(y, 1) - options.firstobs + 1;  % Effective sample size
options.nunits        = size(y, 3);                         % 1 for time series (not panel)

% === PRIOR SETTINGS ===
if prior_on == 1
    options.dummy               = 1;            % Use Minnesota dummy prior
    options.minn_prior_tau      = 3;            % Overall tightness
    options.minn_prior_decay    = 0.5;          % Lag decay
    options.minn_prior_lambda   = 5;            % Sum-of-coefficient prior
    options.minn_prior_mu       = 1;            % Co-persistence
    options.minn_prior_omega    = 0;            % Prior on shock variances
    options.unitroot            = [1; 1; 1; 1; 1; 1];  % Treat all as persistent
else
    options.dummy               = [];
    options.minn_prior_tau      = [];
    options.minn_prior_decay    = [];
    options.minn_prior_lambda   = [];
    options.minn_prior_mu       = [];
    options.minn_prior_omega    = [];
    options.unitroot            = [];
end

% === REGIME SETTINGS ===
if regime_on == 1
    options.regimes = [
        datetime(1973, 1, 1), ...
        datetime(1990, 1, 1)
    ];
else
    options.regimes = [];
end

% === STRUCTURAL IDENTIFICATION SETTINGS ===
options.tvA                   = 0;  % No time-varying A0
options.noLmd                 = 0;  % Estimate shock variances (homoskedastic)
options.A0_restriction        = []; % No A0 prior restrictions
% options.tparam         = 5.703233415698 / 2;  % Scale for A0 prior
% options.tscale         = 1;
options.alpha = 1;
options.K = 1;

% === SIGN RESTRICTIONS (Kilian Style) ===
% Sign restrictions on 3 shocks (each corresponds to a column in IRFs)
options.SignRestrictions = {    
    'y(1,1:1,1) = 0', ...   % zero restriction on the contemporanouos response of TFP (variable 1) to a unit optimism shock
    'y(2,1:1,1) > 0'        % positive sign restriction on the contemporanouos response of Stock Prices (variable 2) to a unit optimism shock
};
info0 = SetupZeroInfo(options.ny, 1, options.SignRestrictions);
options.sign_regime_dependent = (svar_svsr_on == 1);  % Rotational invariance
options.sign_inneriteration = 1000;
options.penalty_offdiagonal_on = (penalty_offdiagonal_on == 1); 

% === DATA STRUCT ===
vardata.y              = y;
vardata.varnames       = varnames;
vardata.filtered_data  = filtered_data;

% === RUN SVAR === 
tStart = tic;
svar_output = svar_run_sign_parallel(vardata, options);

elapsedTime = toc(tStart);
fprintf('Total run time: %.2f seconds (%.2f minutes)\n', elapsedTime, elapsedTime / 60);

%********************************************************
%% Draw impulse responses
%********************************************************
% === Determine IRF output filename based on estimation configuration ===
if regime_on == 0
    % No regime switching: benchmark identification à la Uhlig (2005)
    outname = 'output/SignZero/irfs_uhlig.mat';
else        
    if svar_svsr_on == 1
        % Dawis method: full Q rotation + regime-dependent sign restrictions
        if penalty_offdiagonal_on==1
            outname = 'output/SignZero/irfs_dawis_regime_penalty.mat';
        else
            outname = 'output/SignZero/irfs_dawis_regime.mat';
        end 
    else
        % Carriero et al. (2024): full Q rotation + regime-invariant sign restrictions
        % Using Metropolis-Hastings step
        outname = 'output/SignZero/irfs_carriero_regime.mat';
    end
end

% === Optional suffix modification if Minnesota-style prior is used ===
if prior_on == 1
    [folder, name, ext] = fileparts(outname);
    outname = fullfile(folder, [name, '_minnesota', ext]);
end

save(outname, 'svar_output');
disp(['Saved IRFs to: ', outname]);

system('shutdown /s /t 60'); 