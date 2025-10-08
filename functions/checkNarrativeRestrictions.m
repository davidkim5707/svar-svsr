function isValid = checkNarrativeRestrictions(restrictions, v_narrative)
    % restrictions: cell array of strings like 'v[3,2] > 0.5'
    % v_narrative: numeric matrix

    count = 0;

    for i = 1:length(restrictions)
        r = restrictions{i};

        % Preprocess to replace "v([166],1)" with "v[166,1]"
        r = regexprep(r, 'v\(\[(\d+)\],(\d+)\)', 'v[$1,$2]');
    
        % Extract pattern
        tokens = regexp(r, 'v\[(\d+),(\d+)\]\s*(>|<|>=|<=|==|~=)\s*(-?\d+\.?\d*)', 'tokens');

        if isempty(tokens)
            error(['Invalid restriction format: ', r]);
        end

        tok = tokens{1};
        t     = str2double(tok{1});
        var   = str2double(tok{2});
        op    = tok{3};
        value = str2double(tok{4});

        val = v_narrative(t, var);  % Index from 1 (MATLAB style)
        ok = false;

        switch op
            case '>'
                ok = val > value;
            case '<'
                ok = val < value;
            case '>='
                ok = val >= value;
            case '<='
                ok = val <= value;
            case '=='
                ok = val == value;
            case '~='  % MATLAB uses '~=' for '!='
                ok = val ~= value;
            otherwise
                error(['Unsupported operator: ', op]);
        end

        count = count + ok;
    end

    isValid = (count == length(restrictions));
end
