function w = DrawW(W)
    dim = sum(cellfun(@(A) size(A,1), W));
    w   = zeros(dim,1);
    k=1;
    for j=1:numel(W)
        s = size(W{j},1);
        if s>0
            v = randn(s,1);
            w(k:k+s-1) = v ./ norm(v);
            k = k + s;
        end
    end
end