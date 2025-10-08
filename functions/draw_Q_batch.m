function Qbatch = draw_Q_batch(ny, Ncand)
    Qbatch = zeros(ny, ny, Ncand);
    for k = 1:Ncand
        [Qdraw, R] = qr(randn(ny));
        In    = diag(sign(diag(R)));
        Qdraw = Qdraw  * In;
        Qbatch(:,:,k) = Qdraw;
    end
end