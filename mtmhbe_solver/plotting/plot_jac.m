function p = plot_jac(M, zid, ccolor, msize, xyshift)
    if ~exist('xyshift', 'var')
        xyshift = [];
    end
    if ~exist('msize', 'var')
        msize = [];
    end
    m = M.grid.N;
    n = m;
    A = sparse(M.m, M.n, M.Cpz(:,zid)*M.b_z(zid), n, n);
    pos = get(gca,'position');
    if isempty(msize)
        msize = 0.9*max(4,min(14,round(6*min(pos(3:4))/max(m+1,n+1))));
    end
    [i,j] = find(A);
    if ~isempty(xyshift)
        if numel(xyshift)==2
            i = double(i) + xyshift(1);
            j = double(j) + xyshift(2);
        elseif numel(xyshift)==1
            i = double(i) + xyshift(1);
            j = double(j) + xyshift(1);
        end
    end
    p=scatter3(j, i, A(A~=0), msize, ccolor, 'filled');
end

