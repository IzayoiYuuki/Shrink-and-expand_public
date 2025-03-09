% Test inveres subsapce iteration, SD, LOBPCG, TraceMIN, Chebyshev-Davidson algorithm
% with fix, slope, slopek strategy
% on the test matrices

% (For a quick test, we only test matrix No.1, 2, and 3 here,)
% (which are bcsstm21, rail_5177, and Muu.)
% (If you want to test more, you can set the index of test matrices)
% (in test_'algorithm'_'strategy'.m, respectively.)

clear;
addpath(genpath('./'));

Alg = ["ISI", "SD", "LOBPCG", "TM", "CD"];
Str = ["fix", "slope", "slopek"];

for AlgNo = 1 : 5
    for StrNo = 1 : 3

        disp("Now testing " + Alg(AlgNo) + " with " + Str(StrNo) + " strategy");
        NowTest = "test_" + Alg(AlgNo) + "_" + Str(StrNo);

        try
            exctest(NowTest);
        catch
            disp(Alg(AlgNo) + " with " + Str(StrNo) + " failed!");
        end

    end
end

function [] = exctest(funname)
    eval(funname);
end