function [] = WarmUp(FileName)

disp("Now warming up...");

Maxiter = 3500;
tol = 1e-10;

SEconfig.rule = 'fix';
SEconfig.enlargesteps = 2;
SEconfig.shrinksteps = 10;
warmupiter = 20;
SEconfig.warmuptol = 1e-4;

[A, B] = LoadEigProb(FileName);

% shift the matrix if necessary
    el = eigs(A, 1, 'smallestreal');
if el < 0
    A = A - (1.05*el)*speye(size(A));
end

nev = 100;
nex = ceil(2*nev);

% initial guess
[n, ~] = size(A);
X = randn(n, nex);

% preconditioner
Myprec = @(X) (X);

% SD without shrink
SEconfig.warmupiter = Maxiter;

mySteepestDescent(A, X, Myprec, nev, tol, Maxiter, SEconfig);


% SD with shrink
SEconfig.warmupiter = min(max(round(n/500),SEconfig.shrinksteps), warmupiter);

mySteepestDescent(A, X, Myprec, nev, tol, Maxiter, SEconfig);

hold off;

disp("Warming up finished.");