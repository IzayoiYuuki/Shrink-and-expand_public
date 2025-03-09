function [X] = CholeskyBQR(X, B)
% orthogonalization under B inner product
% output X will satisfying X'BX = I

[~, L] = size(X);

while 1
    [R, flag] = chol(X'*B*X);
    if flag ~= 0
        R = chol(X'*B*X + 1e-10*eye(L));
    end
    X = X/R;
    if norm(X'*B*X - eye(L)) < 1e-11
        break;
    end
end
