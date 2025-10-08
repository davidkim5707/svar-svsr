function W = build_W_from_QZ(Qpart, ZF)
    n = numel(ZF);
    W = cell(1,n);
    for j = 1:n
        Cj = [Qpart(:, 1:max(j-1,0))'; ZF{j}];  % ((j-1)+rows(ZF{j})) x n
        W{j} = null(Cj,'r')';                   % (n - rank(Cj)) x n
    end
end