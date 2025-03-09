function [A, B] = LoadEigProb(filename)
% load A from filename, and load B as identity

load(filename);

A = Problem.A;

[m, n] = size(A);

B = spdiags(ones(n,1), [0], n, n);

