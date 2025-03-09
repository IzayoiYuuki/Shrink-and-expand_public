function [Xr, Dr, iter, res, shrinklist] = myTraceMin(A, B, X, nev, tol, Maxiter, Sconfig, SEconfig)

[n, nex] = size(X);
nowshrink = 0;

normA = normest(A, 0.1);
[X, ~] = qr(X, 0);

shrinklist = zeros(Maxiter, 1);

% Step 0: initial res
AX = A*X;
Ap = X'*AX;
Ap = (Ap' + Ap)/2;
[Z, D] = eig(Ap);
X = X*Z;
AX = AX*Z;
[~, idxk] = sort(diag(D), 'ascend');
X = X(:, idxk); AX = AX(:, idxk); D = D(idxk, idxk); 

% If sparse direct is used, prepare LU decompostion for use
if strcmp(Sconfig.method, 'direct_lu') == 1
    [Sconfig.AL, Sconfig.AU, Sconfig.p1, Sconfig.p2, Sconfig.p3]...
        = prepLU(A);
end

for iter = 1 : Maxiter
    
    % Compute residual
    R = AX - X*D;
    normR = norm(R(:, 1:nev), 'fro');

    % Check convergence
    res_vec = zeros(nex, 1);
    for i = 1 : nex
        res_vec(i) = norm(R(:, i))...
            /(normA*norm(X(:, i)) + norm(X(:, i))*abs(D(i, i)));
    end
    [tmp, idx] = mink(res_vec, nev);
    Xr = X(:, idx);
    Dr = D(idx, idx);
    res(iter) = tmp(end);
    if res(iter) < tol
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

    % Solve PAPd = Pr
    Delta = LagrangeEq(A, B, X, R, Sconfig);

    % Update X
    X = X - Delta;
    [X, ~] = qr(X, 0);

    % Enlarge
    if ifenlarge()
        nowshrink = 0;
        shrinklist(iter) = 1;
        X = [X, Xlog];
        [X, ~] = qr(X, 0);
        [~, nex] = size(X);
    end

    % Solve projection problem
    AX = A*X;
    Ap = X'*AX;
    Ap = (Ap' + Ap)/2;
    [Xp, D] = eig(Ap);
    [~, idx] = sort(abs(diag(D)), 'ascend');
    Xp = Xp(:, idx);
    D = D(idx, idx);
    X = X*Xp;
    AX = AX*Xp;

    % Shrink
    if ifshrink()
        shrinklist(iter) = -1;
        nowshrink = iter;
        nex = nev + 5;
        Xlog = X(:, nex + 1:end);
        X = X(:, 1:nex);
        AX = AX(:, 1:nex);
        D = D(1:nex, 1:nex);
    end
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
        ifelg = nowshrink > 0 && iter - SEconfig.shrinksteps > 0 && shrinklist(iter - SEconfig.shrinksteps) < 0;
    end
end

function [ifshr] = ifshrink()
    ifshr = 0;
    if strcmp(SEconfig.rule, 'slope') || strcmp(SEconfig.rule, 'slopek') || strcmp(SEconfig.rule, 'fix')
        ifshr = (iter >= SEconfig.warmupiter && res(iter) <= SEconfig.warmuptol)...
            && ( ~any(shrinklist) || shrinklist(iter - SEconfig.enlargesteps) > 0);
    end
end

end

function [AL, AU, p1, p2, p3] = prepLU(A)
% If using sparse direct, prepare the LU decompision first
    n = length(A);
    p1 = symamd(A);
    p1 = p1';
    e = speye(n, n);
    e = e(:, p1);
    p3 = e*(1:n)';
    [AL, AU, p2] = lu(A(p1, p1), 'vector');
    p2 = p2';
end

function [Delta] = LagrangeEq(A, B, X, R, Sconfig)
% Solving padel linear system (In sparse direct or iteration methods)
    [n, nex] = size(X);
    Delta = zeros(n, nex);
    if strcmp(Sconfig.method, 'direct') == 1
        Z = A\X;
        Z = Z/(X'*Z);
        Delta = X - Z;
    elseif strcmp(Sconfig.method, 'direct_lu') == 1
        Z = X(Sconfig.p1, :);
        Z = Z(Sconfig.p2, :);
        Z = Sconfig.AL\Z;
        Z = Sconfig.AU\Z;
        Z = Z(Sconfig.p3, :);
        Z = Z/(X'*Z);
        Delta = X - Z;
    else
        R = MKproj(X, R);
        PAPnow = @(z) PAPfunc(A, B, X, z, 0);
        for eqid = 1 : nex
            if strcmp(Sconfig.method, 'pcg') == 1
                % pcg
                [Delta(:, eqid), ~] =...
                    pcg(PAPnow, R(:, eqid), Sconfig.res, Sconfig.iter);
            elseif strcmp(Sconfig.method, 'minres') == 1
                % minres
                [Delta(:, eqid), ~] =...
                    minres(PAPnow, R(:, eqid), Sconfig.res, Sconfig.iter);
            end
        end
    end
end

function [x] = PAPfunc(A, B, X, x, sigma)
% Matrix multiply: x <- P*A*x
% or, precisely, x <- (I - BX*(BX'*BX)^{-1}BX')
%     x = MKproj(BX, x);
    x = A*x - sigma*x;
    x = MKproj(X, x);
end

function [y] = MKproj(X, x)
% Projection on BX: x <- x - BX(BX'BX)^{-1}BX*x
    y = X'*x;
    y = X*y;
    y = x - y;
end
