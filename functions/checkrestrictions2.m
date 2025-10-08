function [d, fsign, diag] = checkrestrictions2(restrictions, y, tol, verbose)
% CHECKRESTRICTIONS2  Evaluate sign & elasticity restrictions without eval.
%   [d, fsign] = checkrestrictions2(restrictions, y)
%   [d, fsign, diag] = checkrestrictions2(restrictions, y, tol, verbose)
%
% Inputs
%   restrictions : cellstr of expressions like:
%       'y(1,1:3,1) < 0'
%       'y(2,2,2) >= 0'
%       'y(1,1:3,2)/y(3,1:3,2) <= 0.0258'
%   y            : (nvar × horizon × nshock) IRFs
%   tol          : numeric tolerance (default 1e-12)
%   verbose      : if true, prints which lines failed to parse/match
%
% Outputs
%   d, fsign     : pass flag and global sign flip used (1 or -1)
%   diag         : diagnostics: pass_vec, kinds, fail indices, counts

if nargin < 3 || isempty(tol), tol = 1e-12; end
if nargin < 4, verbose = false; end

[pass1, diag1] = check_all_restrictions(y, restrictions, tol, verbose);
if pass1, d = 1; fsign = 1; diag = diag1; return; end

[pass2, diag2] = check_all_restrictions(-y, restrictions, tol, verbose);
if pass2, d = 1; fsign = -1; diag = diag2; return; end

d = 0; fsign = 1; diag = diag1;
end


function [all_pass, diag] = check_all_restrictions(y, restrictions, tol, verbose)

n = numel(restrictions);
pass_vec = false(n,1);
kind     = strings(n,1);
fail_ix  = cell(n,1);

sign_needed  = false(n,1);
elastic_seen = false(n,1);
n_sign_total = 0; n_sign_pass = 0;
n_el_total   = 0; n_el_pass   = 0;

% ---- Explicit patterns (no non-capturing groups) ----
% Sign, single t
SIGN_SINGLE = ['^y\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)\s*' ...
               '(>=|<=|>|<)\s*0\s*$'];
% Sign, range t1:t2
SIGN_RANGE  = ['^y\(\s*(\d+)\s*,\s*(\d+)\s*:\s*(\d+)\s*,\s*(\d+)\s*\)\s*' ...
               '(>=|<=|>|<)\s*0\s*$'];

% Elasticity, single t
EL_SINGLE = ['^y\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)\s*/\s*' ...
             'y\(\s*(\d+)\s*,\s*\2\s*,\s*\3\s*\)\s*' ...
             '(<=|>=)\s*' ...
             '([0-9]*\.?[0-9]+([eE][+-]?\d+)?)\s*$'];
% Elasticity, range t1:t2
EL_RANGE  = ['^y\(\s*(\d+)\s*,\s*(\d+)\s*:\s*(\d+)\s*,\s*(\d+)\s*\)\s*/\s*' ...
             'y\(\s*(\d+)\s*,\s*\2\s*:\s*\3\s*,\s*\4\s*\)\s*' ...
             '(<=|>=)\s*' ...
             '([0-9]*\.?[0-9]+([eE][+-]?\d+)?)\s*$'];

for k = 1:n
    restr = strtrim(restrictions{k});

    % ---------- SIGN: single ----------
    tok = regexp(restr, SIGN_SINGLE, 'tokens');
    if ~isempty(tok)
        kind(k) = "sign"; sign_needed(k) = true; n_sign_total = n_sign_total + 1;
        t = tok{1};
        i = str2double(t{1}); t1 = str2double(t{2}); t2 = t1; s = str2double(t{3}); op = t{4};
        vals = reshape(y(i, t1:t2, s), 1, []);
        ok = cmp_sign(vals, op, tol);
        if ok, pass_vec(k)=true; n_sign_pass = n_sign_pass + 1; else, fail_ix{k}=t1:t2; end
        continue
    end

    % ---------- SIGN: range ----------
    tok = regexp(restr, SIGN_RANGE, 'tokens');
    if ~isempty(tok)
        kind(k) = "sign"; sign_needed(k) = true; n_sign_total = n_sign_total + 1;
        t = tok{1};
        i = str2double(t{1}); t1 = str2double(t{2}); t2 = str2double(t{3}); s = str2double(t{4}); op = t{5};
        vals = reshape(y(i, t1:t2, s), 1, []);
        ok = cmp_sign(vals, op, tol);
        if ok, pass_vec(k)=true; n_sign_pass = n_sign_pass + 1; else, fail_ix{k}=t1:t2; end
        continue
    end

    % ---------- ELASTICITY: single ----------
    tok = regexp(restr, EL_SINGLE, 'tokens');
    if ~isempty(tok)
        kind(k) = "elasticity"; elastic_seen(k) = true; n_el_total = n_el_total + 1;
        t = tok{1};
        i = str2double(t{1}); tt = str2double(t{2}); s = str2double(t{3});
        j = str2double(t{4}); op = t{5}; c = str2double(t{6});
        num = reshape(y(i, tt, s), 1, []);
        den = reshape(y(j, tt, s), 1, []);
        ok = cmp_elastic(num, den, op, c, tol);
        if ok, pass_vec(k)=true; n_el_pass = n_el_pass + 1; else, fail_ix{k}=tt; end
        continue
    end

    % ---------- ELASTICITY: range ----------
    tok = regexp(restr, EL_RANGE, 'tokens');
    if ~isempty(tok)
        kind(k) = "elasticity"; elastic_seen(k) = true; n_el_total = n_el_total + 1;
        t = tok{1};
        i = str2double(t{1}); t1 = str2double(t{2}); t2 = str2double(t{3}); s = str2double(t{4});
        j = str2double(t{5}); op = t{6}; c = str2double(t{7});
        num = reshape(y(i, t1:t2, s), 1, []);
        den = reshape(y(j, t1:t2, s), 1, []);
        ok = cmp_elastic(num, den, op, c, tol);
        if ok, pass_vec(k)=true; n_el_pass = n_el_pass + 1; else, fail_ix{k}=t1:t2; end
        continue
    end

    % ---------- Unknown pattern ----------
    kind(k) = "unknown";
    pass_vec(k) = false;
    fail_ix{k} = [];
    if verbose
        fprintf('[parse] Unknown or unsupported restriction: %s\n', restr);
    end
end

% Final decision
if ~any(elastic_seen)
    all_pass = (n_sign_pass == n_sign_total) && all(pass_vec | ~sign_needed);
else
    all_pass = (n_sign_pass == n_sign_total) && (n_el_pass == n_el_total);
end

diag.pass_vec       = pass_vec;
diag.kind           = kind;
diag.fail_ix        = fail_ix;
diag.n_sign_total   = n_sign_total;
diag.n_sign_pass    = n_sign_pass;
diag.n_el_total     = n_el_total;
diag.n_el_pass      = n_el_pass;

end


% ===== helpers =====
function ok = cmp_sign(vals, op, tol)
switch op
    case '>' , ok = all(vals > 0 + tol);
    case '>=' , ok = all(vals >= 0 - tol);
    case '<' , ok = all(vals < 0 - tol);
    case '<=' , ok = all(vals <= 0 + tol);
    otherwise, ok = false;
end
end

function ok = cmp_elastic(num, den, op, c, tol)
safe = all(abs(den) > tol);
if ~safe, ok = false; return; end
ratio = num ./ den;
switch op
    case '<=' , ok = all(isfinite(ratio) & ratio <= c + tol);
    case '>=' , ok = all(isfinite(ratio) & ratio >= c - tol);
    otherwise, ok = false;
end
end
