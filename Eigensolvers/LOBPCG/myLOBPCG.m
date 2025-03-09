function [X, Theta, iter, res, shrinklist] = myLOBPCG(A, B, X,...
    Prec, nev, tol, Maxiter, SEconfig)
% LOBPCG algorithm with soft-locking and HL trick
% computing smallest eigenpairs by magnitude
% for reference, read doi: 10.1137/17M1129830
%           Parameters:
% A, B:    eigenvalue problem - Ax = lambda*Bx
% X:       initial guess
% Prec:    preconditioner
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

% Initialize some data
nowshrink = 0;
Xlog = [];

[~, nex] = size(X);

normA = normest(A, 1e-2);
normB = normest(B, 1e-2);

X = CholeskyBQR(X, B);

[Z, Theta] = eig(X'*A*X);
[~, idx] = sort(diag(Theta), 'ascend');
Z = Z(:, idx); Theta = Theta(idx, idx);

X = X*Z;
P = [];

% Loop
for iter = 1 : Maxiter

    % save shrink log
    shrinklist(iter) = 0;
    
    % residual
    R = A*X - B*X*Theta;
        
    % check convergence
    reslist = zeros(1, nex);
    nconv = 0; stopconv = 0;
    nlock = 0; stoplock = 0;
    for xnow = 1 : nev
       resnow = norm(R(:, xnow))/...
                ((normA + normB*abs(Theta(xnow, xnow)))*norm(X(:, xnow)));
       reslist(xnow) = resnow;
       if stopconv ~= 1 && resnow < tol
           nconv = nconv + 1;
           if resnow < tol/100 && stoplock ~= 1
               nlock = nlock + 1;
           else
               stoplock = 1;
           end
       else
           stopconv = 1;
       end
    end
    
    res(iter) = max(reslist);

    % sucessful exit
    if nconv >= nev
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
    
    % deflation
    if iter == 1
       R = R(:, end - nex + nlock + 1:end);
    else
       [~, np] = size(P);
       vecwant = nex - nlock;
       vecwant = min(np, vecwant);
       R = R(:, end - vecwant + 1:end);
       P = P(:, end - vecwant + 1:end);
    end
    
    % precondition
    R = Prec(R);
    
    % Ortho
    XP = [X, P];
    R = R - XP*(XP'*(B*R));
    R = CholeskyBQR(R, B);
    R = R - XP*(XP'*(B*R));
    R = CholeskyBQR(R, B);
    
    S = [X, P, R];
    
    % enlarge the subsapce
    if ifenlarge()
       % save the shrink log
       shrinklist(iter) = +1;
       % mark that enlarge has been done
       nowshrink = 0;
       [~, Ln] = size(Xlog);
       nex = nex + Ln;
    
       Selg = [Xlog, Plog];
       Selg = Selg - S*(S'*(B*Selg));
       Selg = CholeskyBQR(Selg, B);
       
       X = [X, Selg(:, 1:end/2)];
       P = [P, Selg(:, end/2 + 1:end)];
       S = [X, P, R];
    end
    
    % Rayleigh-Ritz
    Ap = S'*A*S;
    Ap = (Ap + Ap')/2;
    [Z, Theta] = eig(Ap);
    [~, idx] = sort(diag(Theta), 'ascend');
    Z = Z(:, idx);
    Theta = Theta(idx, idx); Theta = Theta(1:nex, 1:nex);
    
    % HL trick
    if nex <= length(Z)/2
       Cp = Z(1:nex, nex + 1:end)'*Z(1:nex, nlock + 1:nex);
    else
       Cp = Z(nex + 1:end, nex + 1:end)'*Z(nex + 1:end, nlock + 1:nex);
    end
    [Cp, ~] = qr(Cp, 0);
    Cp = Z(:, nex + 1:end)*Cp;
    Cx = Z(:, 1:nex);
    
    % update X P
    X = S*Cx;
    P = S*Cp;
    
    % shrink the subspace
    if ifshrink()
       % save the shrink log
       shrinklist(iter) = -1;
       % mark that shrink has been down
       nowshrink = iter;
       % shrink to nev + 5
       Sn = nev + 5 - nex;
       nex = nex + Sn;
       % save X and P
       Xlog = X(:, end + Sn + 1:end);
       Plog = P(:, end + Sn + 1:end);
       X = X(:, 1:end + Sn);
       P = P(:, 1:end + Sn);
       Theta = Theta(1:nex, 1:nex);
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

