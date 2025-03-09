function [V, D, iter, res, shrinklist] = mySubspaceIteration(A, X, nev, tol, Maxiter, SEconfig)
% a classical subspace iteration algorithm
% computing largest eigenpairs by magnitude
% search space is updated by A*X in each iteration
%           Parameters:
% A:       eigenvalue problem - Ax = lambda*x
% X:       initial guess
% nev:     number of eigenvalue to be computed
% tol:     convergence tolerance of residual norm(Ax-lambda*Bx)/(norm(x)*(norm(A)+abs(lambda)*norm(B)))
% Maxiter: max iteration
% SEconfig: setting of shrink-and-expand technique:
%           SEconfig.rule:         SE strategy
%           SEconfig.enlargesteps: step of expand --> shrink
%           SEconfig.shrinksteps:  step of shrink --> expand (if fix)
%           SEconfig.enlargetol:   tolerance of employing expand (if slope or slopek)
%           SEconfig.slopestep:    number of iteration for taking average (if slopek)
%           SEconfig.warmupiter:   minimum iteration before employing SE
%           SEconfig.warmuptol:    maximum residual before employing SE

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