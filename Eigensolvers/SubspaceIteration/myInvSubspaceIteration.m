function [V, D, iter, res, shrinklist] = myInvSubspaceIteration(A, X, nev, tol, Maxiter, SEconfig)
% a classical inverse subspace iteration algorithm
% computing smallest eigenpairs by magnitude
% search space is updated by A^{-1}X in each iteration
%           Parameters:
% A:       eigenvalue problem - Ax = lambda*x
% X:       initial guess
% nev:     number of eigenvalue to be computed
% tol:     convergence tolerance of residual norm(Ax-lambda*Bx)/(norm(x)*(norm(A)+abs(lambda)*norm(B)))
% Maxiter: max iteration
% SEconfig: setting of shrink-and-expand technique:
%           SEconfig.rule:         SE strategy
%           SEconfig.expandsteps: step of expand --> shrink
%           SEconfig.shrinksteps:  step of shrink --> expand (if fix)
%           SEconfig.expandtol:   tolerance of employing expand (if slope or slopek)
%           SEconfig.slopestep:    number of iteration for taking average (if slopek)
%           SEconfig.warmupiter:   minimum iteration before employing SE
%           SEconfig.warmuptol:    maximum residual before employing SE

normA = normest(A, 1e-3);
shrinklist = zeros(Maxiter, 1);
nowshrink = 0;
[X, ~] = qr(X, 0);

% Prepare LU decomposition
facA = decomposition(A);

% Step 0: initial res
AX = A*X;
Ap = X'*AX;
Ap = (Ap' + Ap)/2;
[Z, D] = eig(Ap);
X = X*Z;
AX = AX*Z;
[~, idxk] = sort(diag(D), 'ascend');
X = X(:, idxk); AX = AX(:, idxk); D = D(idxk, idxk); 

% Loop: Subspace Iteration
for iter = 1 : Maxiter

    % check convergence
    V = X(:, 1:nev); AV = AX(:, 1:nev); D = D(1:nev, 1:nev); 
    R = AV - V*D;
    res_vec = zeros(nev, 1);
    for i = 1 : nev
        res_vec(i) = norm(R(:, i))/(normA + abs(D(i)));
    end
    res(iter) = max(res_vec);
    if (res(iter) <= tol)
        break;
    end

    % To find the steepest decreacion after last shrink
    if strcmp(SEconfig.rule, 'slope') && nowshrink > 0
        nowD = log10(res(nowshrink)) - log10(res(nowshrink + 1));
        newD = log10(res(iter - 1)) - log10(res(iter));
        if newD > nowD
            nowshrink = iter - 1;
        end
    elseif strcmp(SEconfig.rule, 'slopek') && nowshrink > 0 && iter - SEconfig.slopestep > nowshrink
        nowD = log10(res(nowshrink)) - log10(res(nowshrink + SEconfig.slopestep));
        newD = log10(res(iter - SEconfig.slopestep)) - log10(res(iter));
        if newD > nowD
            nowshrink = iter - SEconfig.slopestep;
        end
    end

    % Inverse iteration
    X = facA\X;

    if ifexpand()
        X = [X, Xlog];
        shrinklist(iter) = 1;
        nowshrink = 0;
    end

    % R-R projection
    [X, ~] = qr(X, 0);
    AX = A*X;
    Ap = X'*AX;
    Ap = (Ap' + Ap)/2;
    [Z, D] = eig(Ap);
    X = X*Z;
    AX = AX*Z;

    % choose the best eigen pairs
    [~, idxk] = sort(diag(D), 'ascend');
    X = X(:, idxk); AX = AX(:, idxk); D = D(idxk, idxk);
    
    if ifshrink()
        Xlog = X(:, nev + 6:end);
        X = X(:, 1:nev + 5);
        AX = AX(:, 1:nev + 5);
        shrinklist(iter) = -1;
        nowshrink = iter;
    end

end

function [ifelg] = ifexpand()
    ifelg = 0;
    if strcmp(SEconfig.rule, 'slope')
        ifelg = nowshrink > 0 && (log10(res(nowshrink)) - log10(res(nowshrink + 1)))*(iter - nowshrink)/...
                (log10(res(nowshrink)) - log10(res(iter))) > SEconfig.expandtol;
    elseif strcmp(SEconfig.rule, 'slopek')
        ifelg = nowshrink > 0 && iter - nowshrink >= SEconfig.slopestep...
                && (log10(res(nowshrink)) - log10(res(nowshrink + SEconfig.slopestep)))*(iter - nowshrink)/...
                ((log10(res(nowshrink)) - log10(res(iter)))*SEconfig.slopestep) > SEconfig.expandtol;
    elseif strcmp(SEconfig.rule, 'fix')
        ifelg = nowshrink > 0 && iter - SEconfig.shrinksteps > 0 && shrinklist(iter - SEconfig.shrinksteps) < 0;
    end
end

function [ifshr] = ifshrink()
    ifshr = 0;
    if strcmp(SEconfig.rule, 'slope') || strcmp(SEconfig.rule, 'slopek') || strcmp(SEconfig.rule, 'fix')
        ifshr = (iter >= SEconfig.warmupiter && res(iter) <= SEconfig.warmuptol)...
            && ( ~any(shrinklist) || shrinklist(iter - SEconfig.expandsteps) > 0);
    end
end

end
    
