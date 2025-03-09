function [V, D, iter, res, logs] = myChebyshevDavidson(A, X, nev, tol, Maxiter, CDconfig, SEconfig)

% problem information
[n, nex] = size(X);
normA = normest(A, 1e-1);

% initialize parameters
logs = struct();
logs.restart = [];
logs.shrinklist = zeros(Maxiter, 1);
V = [];
D = [];
nconv = 0;
nowshrink = 0;
nex_expand = nex;
nex_shrink = ceil(nex/2);

% parameter for Cheby filter
lowb = CDconfig.lowb;
upb = CDconfig.upb;
polyorder = CDconfig.polyorder;
submax = CDconfig.submax;
newsub = CDconfig.newsub;

% iteration 1
X = Chebyshev_filter(polyorder, lowb, upb, X, A);

% initialize X
[X, ~] = qr(X, 0);

% Loop
for iter = 1 : Maxiter

    % R-R projection
    AX = A*X;
    Ap = X'*AX;
    Ap = (Ap' + Ap)/2;
    [Z, D] = eig(Ap);
    X = X*Z;
    AX = AX*Z;

    % choose the best eigen pairs
    [~, idxk] = sort(diag(D), 'ascend');
    X = X(:, idxk); AX = AX(:, idxk); D = D(idxk, idxk);

    % residual
    [~, ksub] = size(X); 
    R = AX - X*D;
    res_vec = zeros(min(nev, ksub), 1);
    for i = 1 : min(nev, ksub)
        res_vec(i) = norm(R(:, i))/(normA + abs(D(i)));
    end

    % check convergence
    if ksub < nev
        res(iter) = 1;
    else
        res(iter) = max(res_vec);
        if (res(iter) <= tol)
            V = X(:, 1:nev); D = D(1:nev, 1:nev);
            break;
        end
    end
    [~, idx] = sort(res_vec, 'ascend');
    logs.res_now(iter) = res_vec(idx(nconv + 1));
    logs.nconv(iter) = nconv;
    X(:, 1:length(idx)) = X(:, idx);
    nconv = sum(res_vec < tol);

    % Restart
    if ksub > submax
        logs.restart = [logs.restart, iter];
        % all converged vectors + nex other vectors
        X = X(:, 1:max(nconv + newsub, nev));
    end

    if ifenlarge()
        nex = nex_expand;
        logs.shrinklist(iter) = 1;
        nowshrink = 0;
    elseif ifshrink()
        nex = nex_shrink;
        logs.shrinklist(iter) = -1;
        nowshrink = iter;
    end

    % Chebyshev filter
    Xtmp = X(:, (nconv + 1):(nconv + nex));
    Xtmp = Chebyshev_filter(polyorder, lowb, upb, Xtmp, A);

    % orthogonalization
    Xtmp = Xtmp - X*(X'*Xtmp);
    [Xtmp, ~] = qr(Xtmp, 0);
    Xtmp = Xtmp - X*(X'*Xtmp);
    [Xtmp, ~] = qr(Xtmp, 0);
    X = [X, Xtmp];

    % X = [X, Xtmp];
    % [X, ~] = qr(X, 0);
    
end
    
function [ifelg] = ifenlarge()
    ifelg = 0;
    if strcmp(SEconfig.rule, 'slope')
        ifelg = nowshrink > 0 && (log10(res(nowshrink)) - log10(res(nowshrink + 1)))*(iter - nowshrink)/...
                (log10(res(nowshrink)) - log10(res(iter))) > SEconfig.enlargetol;
    elseif strcmp(SEconfig.rule, 'slopek')
        ifelg = nowshrink > 0 && iter - nowshrink >= SEconfig.slopestep...
                && (log10(res(nowshrink)) - log10(res(nowshrink + SEconfig.slopestep)))*(iter - nowshrink)/...
                ((log10(res(nowshrink)) - log10(res(iter)))*SEconfig.slopestep) > SEconfig.enlargetol;
    elseif strcmp(SEconfig.rule, 'fix')
        ifelg = nowshrink > 0 && iter - SEconfig.shrinksteps > 0 && logs.shrinklist(iter - SEconfig.shrinksteps) < 0;
    end
end

function [ifshr] = ifshrink()
    ifshr = 0;
    if strcmp(SEconfig.rule, 'slope') || strcmp(SEconfig.rule, 'slopek') || strcmp(SEconfig.rule, 'fix')
        ifshr = (iter >= SEconfig.warmupiter && res(iter) <= SEconfig.warmuptol)...
            && ( ~any(logs.shrinklist) || logs.shrinklist(iter - SEconfig.enlargesteps) > 0);
    end
end

end

function [Y] = Chebyshev_filter(m, a, b, X, A)
    % Filter X by an degree Chebyshev polynomial
    e = (b - a)/2;
    c = (b + a)/2;
    Y = (A*X - c*X)/e;
    
    for i = 2 : m
        Ynew = 2*(A*X - c*X)/e - X;
        X = Y;
        Y = Ynew;
    end
end