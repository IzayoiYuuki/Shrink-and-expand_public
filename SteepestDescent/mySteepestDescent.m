function [X, Theta, iter, res, shrinklist] = mySteepestDescent(A, X, precond, nev, tol, maxiter, SEconfig)
% A should be symmetric
nex = size(X, 2);
shrinklist = zeros(maxiter, 1);
nowshrink = 0;
Xlog = [];

normA = normest(A, 1e-2);

[X, ~] = qr(X, 0);

[C, Theta] = eig(X'*A*X);
[~, idx] = sort(abs(diag(Theta)), 'ascend');
C = C(:, idx); Theta = Theta(idx, idx);

X = X*C;
endwarmup = 0;

for iter = 1 : maxiter
    % save shrink log
    shrinklist(iter) = 0;

    % residual
    R = A*X - X*Theta;

    % check convergence
    reslist = zeros(1, nex);
    for xnow = 1 : nex
       resnow = norm(R(:, xnow))/...
                ((normA + abs(Theta(xnow, xnow)))*norm(X(:, xnow)));
       reslist(xnow) = resnow;
    end
    reslist = mink(reslist, nev);
    res(iter) = reslist(end);

    % sucessful exit
    if res(iter) <= tol
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


    % enlarge the subsapce
    if ifenlarge()
        % save the shrink log
        shrinklist(iter) = +1;
        % mark that enlarge has been done
        nowshrink = 0;
        [~, Ln] = size(Xlog);
        nex = nex + Ln;

        X = [X, Xlog];
    end

    % Projection subspace
    [S, ~] = qr([X, R], 0);

    % Rayleigh-Ritz
    Ap = S'*A*S;
    Ap = (Ap + Ap')/2.0;
    [C, Theta] = eig(Ap);
    [~, idx] = sort(diag(Theta), 'ascend');
    C = C(:, idx); Theta = Theta(idx, idx);
    X = S*C(:, 1:nex); Theta = Theta(1:nex, 1:nex);

    % shrink the subspace
    if ifshrink()
        % save the shrink log
        shrinklist(iter) = -1;
        % mark that shrink has been down
        nowshrink = iter;
        % shrink to nev + 5
        Sn = nev + 5 - nex;
        nex = nex + Sn;
        % save X and R
        Xlog = X(:, end + Sn + 1:end);
        X = X(:, 1:end + Sn);
        Theta = Theta(1:end + Sn, 1:end + Sn);
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
    if iter >= SEconfig.warmupiter && res(iter) <= SEconfig.warmuptol
        endwarmup = 1;
    end
    if strcmp(SEconfig.rule, 'fix') || strcmp(SEconfig.rule, 'slope') || strcmp(SEconfig.rule, 'slopek')
        ifshr = (endwarmup == 1) && ( ~any(shrinklist) || shrinklist(iter - SEconfig.enlargesteps) > 0);
    end
end

end