function [V, D, iter, res, shrinklist] = mySubspaceIteration(A, X, nev, tol, Maxiter, SEconfig)

normA = normest(A, 1e-3);
shrinklist = zeros(Maxiter, 1);
nowshrink = 0;
[X, ~] = qr(X, 0);

% Step 0: initial res
[V, D] = eig(X'*A*X);
X = X*V;
[~, idxk] = maxk(diag(D), nev);
res(1) = norm(A*X(:, idxk) - X(:, idxk)*D(idxk, idxk))/normA;

% Loop: Subspace Iteration
for iter = 1 : Maxiter
    [X, ~] = qr(A*X, 0);
    [V, D] = eig(X'*A*X);
    X = X*V;

    % choose the best eigen pairs
%     [~, idxk] = maxk(diag(D), nev);
    [~, idxk] = sort(diag(D), 'descend');
    X = X(:, idxk);
    D = D(idxk, idxk);
    
    % save temporary result
    V = X(:, 1:nev);
    D = D(1:nev, 1:nev);

    % check convergence
    R = A*V - V*D;
    res_vec = zeros(nev, 1);
    for i = 1 : nev
        res_vec(i) = norm(R(:, i))/(normA + abs(D(i)));
    end
    res(iter + 1) = max(res_vec);
    if (res(iter + 1) <= tol)
        break;
    end

    if ifshrink()
        Xlog = X(:, nev + 1:end);
        X = X(:, 1:nev);
        shrinklist(iter + 1) = -1;
        nowshrink = 1;
    end

    if ifenlarge()
        X = [X, Xlog];
        [X, ~] = qr(X, 0);
        shrinklist(iter + 1) = 1;
        nowshrink = 0;
    end

end

function [ifelg] = ifenlarge()
    ifelg = nowshrink > 0 && iter - SEconfig.shrinksteps > 0 && shrinklist(iter - SEconfig.shrinksteps) < 0;
end

function [ifshr] = ifshrink()
    ifshr = 0;
    ifshr = iter >= SEconfig.warmupiter...
                && ( ~any(shrinklist) || (iter - SEconfig.shrinksteps > 0 && shrinklist(iter - SEconfig.enlargesteps) > 0));
end

end