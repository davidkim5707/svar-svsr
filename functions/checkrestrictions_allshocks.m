function [d, accepted_shock, fsign] = checkrestrictions_allshocks(SignRestrictions, irfs)
%--------------------------------------------------------------------------
% checkrestrictions_allshocks
%   Evaluates sign restrictions for each structural shock in a VAR.
%   Returns success if at least one shock (possibly flipped) satisfies all restrictions.
%
% Inputs:
%   - restrictions     : cell array of strings. Each string should be a logical
%                        expression involving IRFs (evaluated on variable `yj`).
%                        Example: {'all(yj(2,1:4) < 0)', 'all(yj(4,1:2) > 0)'}
%
%   - y                : IRF array of size [nvar x nhorizon x nshock],
%                        where y(i,t,j) is the response of variable i at horizon t to shock j
%
% Outputs:
%   - d                : scalar = 1 if at least one shock satisfies the restrictions
%                                 (either as is or after sign flip), else 0
%   - accepted_shock   : index of the first shock that satisfies the restrictions (NaN if none)
%   - fsign            : +1 if restrictions were satisfied without sign flip,
%                        -1 if sign flip was necessary (NaN if no shock accepted)
%--------------------------------------------------------------------------

nshock = size(irfs, 3);
d = 0;                        % default to reject
accepted_shock = NaN;         % default to no shock accepted
fsign = NaN;                  % default to no sign assigned

% Loop over shocks
for j = 1:nshock
    y = irfs(:,:,j);            % extract IRFs for shock j
    r = SignRestrictions;         % apply same set of restrictions to all shocks

    % --- Check if the IRFs satisfy restrictions without flipping
    count = 0;
    for ii = 1:length(r)
        tmp = eval(r{ii});    % evaluate logical string (e.g., 'all(yj(2,1:4) < 0)')
        count = count + min(tmp);  % count = number of satisfied restrictions
    end
    if count == length(r)
        d = 1;
        accepted_shock = j;
        fsign = +1;
        return;               % success without flipping
    end

    % --- Check if the IRFs satisfy restrictions after flipping sign
    y = -irfs(:,:,j);           % flip sign of all responses for shock j
    count_flip = 0;
    for ii = 1:length(r)
        tmp = eval(r{ii});
        count_flip = count_flip + min(tmp);
    end
    if count_flip == length(r)
        d = 1;
        accepted_shock = j;
        fsign = -1;
        return;               % success with sign flip
    end
end

% If no shock satisfies the restriction set (even after flipping), d = 0
end
