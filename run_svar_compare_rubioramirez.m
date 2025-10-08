%********************************************************
%% Data loading and setting
%********************************************************
clc;        % Clear command window
clear;      % Clear all variables from workspace
close all;  % Close all figures

% Add the "functions" folder to MATLAB path
addpath(genpath('functions'));

% Load data from CSV file
data = load('data/Uhlig_Data_Updated.mat');

% Extract variable names from column 2 to the last column
varnames = data.varNames;

% % Convert date column to datetime format
filtered_data = data;
filtered_data.dates = datetime(data.dates, 'ConvertFrom', 'datenum');
% 
% % Filter data to keep only observations within the date range
% filtered_data = data(find(data.dates == datetime('1965-01-01')):find(data.dates == datetime('2007-12-01')),:);
% filtered_data.dates = datetime(filtered_data.dates, 'InputFormat', 'yyyy-MM-dd');

% Replace negative values with 1 to avoid complex values when taking log
%filtered_data.NONBORRES(filtered_data.NONBORRES <= 0) = 1;

% Drop the 'date' column and define y
y = data.data;

% Convert to array if needed
%y = table2array(y);

% Export to CSV
%writetable(filtered_data(:, [1,2,3,4,5,6,7]), '../../data/filtered_data.csv');

%********************************************************
%% Define Options Struct for SVAR
%********************************************************

% === USER SETTINGS ===
% --------------------------------------------------------------
% Uhlig (2005):             regime_on = 0, prior_on = 0, svar_svsr_on = 0
% Carriero et al. (2024):   regime_on = 1, prior_on = 0, svar_svsr_on = 0
% Dawis' method:            regime_on = 1, prior_on = 0, svar_svsr_on = 1
% --------------------------------------------------------------
regime_on     = 1;    % 1 = enable regime-switching
prior_on      = 0;    % 1 = use Minnesota prior
svar_svsr_on  = 1;    % 1 = Dawis (rotational invariance)
penalty_offdiagonal_on = 1; % 1 = Use penalty function to check sign-restrictions

% === MODEL STRUCTURE ===
options.ny            = size(y, 2);       % Number of observable variables
options.lags          = 12;               % Number of lags
options.irf_horizon   = 60;               % Horizon for IRFs
options.constant      = 1;                % Include constant in VAR
options.timetrend     = 0;                % Include time trend
options.nexogenous    = 0;                % No exogenous variables
options.exogenous     = [];               % Empty exogenous matrix
options.non_explosive_ = 0;               % Do not enforce non-explosiveness

% === SAMPLING SETTINGS ===
options.firstobs      = options.lags + 1;                     % First observation for estimation
options.presample     = 0;                                    % Number of presample periods
options.ndraw         = 3000;                                 % Total recorded draws
options.nsep          = 1;                                    % Draws per stored draw
options.nburn         = 1000;                                 % Burn-in
options.max_compute   = 1;                                    % Optimization algorithm (1 = csminwel)
options.nit           = 50;                                   % Max iterations
options.hsnscale      = 0.01;                                  % Hessian scaling factor

% === DATA DIMENSIONS ===
options.nobs          = size(y, 1) - options.firstobs + 1;    % Usable observations
options.nunits        = size(y, 3);                           % Number of units (1 if time series)

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
        datetime(1979, 10, 1), ...  % Volcker shock: monetarist experiment begins
        datetime(1983, 1, 1),  ...  % End of monetarist period
        datetime(1990, 1, 1)
    ];
else
    options.regimes = [];
end

% === STRUCTURAL IDENTIFICATION SETTINGS ===
options.tvA                   = 0;  % No time-varying A0
options.noLmd                 = 0;  % Estimate shock variances (homoskedastic)
options.A0_restriction = []; % options.A0_restriction = triu(true(options.ny));  % upper triangular: diagonal + upper part = true
options.tparam         = 5.703233415698 / 2;  % Scale for A0 prior
options.tscale         = 1;
% options.alpha = 1;
% options.K = 1;

% === SIGN RESTRICTIONS ===
options.sign_horizon = 6;
options.SignRestrictions = {
    sprintf('y(2,1:%d,1) < 0', options.sign_horizon), ...
    sprintf('y(3,1:%d,1) < 0', options.sign_horizon), ...
    sprintf('y(5,1:%d,1) < 0', options.sign_horizon), ...  % Interest rate ↑
    sprintf('y(6,1:%d,1) > 0', options.sign_horizon)
};
options.sign_regime_dependent = (svar_svsr_on == 1);  % Rotational invariance
options.sign_inneriteration = 1000;
options.penalty_offdiagonal_on = (penalty_offdiagonal_on == 1); 

% === NARRATIVE RESTRICTIONS ===
options.NarrativeRestrictions.time_index = 166;     % time (observation index)
options.NarrativeRestrictions.shock_index = 1;      % shock (column in structural shocks)
options.NarrativeRestrictions.sign = 1;   
options.NarrativeRestrictions.variable_index = 6;

% === DATA STRUCT ===
vardata.y              = y;
vardata.varnames       = varnames;
vardata.filtered_data  = filtered_data;

% === RUN SVAR ===
svar_output = svar_run_sign_parallel(vardata, options);

% %********************************************************
% %% Draw impulse responses
% %********************************************************
% === Determine IRF output filename based on estimation configuration ===

if regime_on == 0
    % === No regime-switching: benchmark SVAR à la Uhlig (2005) ===
    outname = 'output/irfs_narrative_rubioramirez.mat';

else
    % === Regime-switching model activated ===

    if svar_svsr_on == 1
        % === Dawis method: regime-dependent sign restrictions on Q ===
        if penalty_offdiagonal_on == 1
            % Penalty term added for off-diagonal elements of Λ
            outname = 'output/irfs_narrative_dawis_regime_penalty.mat';
        else
            % Standard Dawis regime-switching model without penalty
            outname = 'output/irfs_narrative_dawis_regime.mat';
        end
    else
        % === Carriero et al. (2024): regime-invariant sign restrictions ===

        outname = 'output/irfs_narrative_carriero_regime.mat';
    end
end

% === Optional suffix modification if Minnesota-style prior is used ===
if prior_on == 1
    [folder, name, ext] = fileparts(outname);
    outname = fullfile(folder, [name, '_minnesota', ext]);
end

% % === If the file already exists, append _ver2 before the extension ===
% if exist(outname, 'file')
%     [folder, name, ext] = fileparts(outname);
%     outname = fullfile(folder, [name, '_ver2', ext]);
% end

save(outname, 'svar_output');
disp(['Saved IRFs to: ', outname]);

system('shutdown -s -t 60'); 