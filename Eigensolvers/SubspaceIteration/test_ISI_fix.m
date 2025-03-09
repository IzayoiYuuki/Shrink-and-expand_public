clear;
hold off;

% Add root path
ScriptPath = fileparts(mfilename("fullpath"));
RootPath = fullfile(ScriptPath, "../../");
addpath(genpath(RootPath));

% Add matrix path & save path
MatNames = ["bcsstm21";
            "rail_5177";
            "Muu";
            "fv1";
            "shuttle_eddy";
            "barth5";
            "Si5H12";
            "mario001";
            "c-65";
            "Andrews";
            "Ga3As3H12";
            "Ga10As10H30"];
FileNames = fullfile(RootPath, "Matrices", MatNames + ".mat");
FigureOut = fullfile(RootPath, "Figure", MatNames + "_ISI_fix.pdf");
DataOut = fullfile(RootPath, "Figure", "Data_ISI_fix.txt");

% Problem setting
Nevs = [100; 100; 100; 100; 104; 156; 199; 384; 481; 500; 500; 500];
Maxiter = 1200;
tol = 1e-10;

% solver setting
SEconfig.rule = 'fix';
SEconfig.enlargesteps = 2;
SEconfig.shrinksteps = 10;
warmupiter = 5;
SEconfig.warmuptol = 1e-4;

% Warm up
WarmUp(FileNames(2));

dlmwrite(DataOut, date, '-append', 'delimiter', '', 'precision', 4);

for fileNo = 1 : 3

    disp(MatNames(fileNo));
    dlmwrite(DataOut, fileNo, '-append', 'delimiter', '', 'precision', 4);

    [A, B] = LoadEigProb(FileNames(fileNo));

    % shift the matrix if necessary
    el = eigs(A, 1, 'smallestreal');
    if el < 0
        A = A - (1.05*el)*speye(size(A));
    end

    nev = Nevs(fileNo);
    nex = ceil(2*nev);

    % initial guess
    rng(0);
    [n, ~] = size(A);
    X = randn(n, nex);

    % preconditioner
    Myprec = @(X) (X);

    % ISI without shrink
    SEconfig.warmupiter = Maxiter;
    tic;
    [~, ~, iter, res, ~] = myInvSubspaceIteration(A, X, nev, tol, Maxiter, SEconfig);
    timeL(fileNo, 1) = toc;
    iterL(fileNo, 1) = iter;
    resL{fileNo}(1, 1:length(res)) = res;

    semilogy((1:iter), res, '-*', 'linewidth', 2);
    hold on;

    % ISI with shrink
    SEconfig.warmupiter = warmupiter;
    tic;
    [~, ~, iter, res, SEconfig.shrinklist] = myInvSubspaceIteration(A, X, nev, tol, Maxiter, SEconfig);
    timeL(fileNo, 2) = toc;
    iterL(fileNo, 2) = iter;
    resL{fileNo}(2, 1:length(res)) = res;
    shrinklistL(fileNo, 1:length(SEconfig.shrinklist)) = SEconfig.shrinklist;

    semilogy((1:iter), res, '-o', 'linewidth', 2); hold on;
    SEconfig.shrinklist = SEconfig.shrinklist(1:iter);
    plot(find(SEconfig.shrinklist < 0), res(SEconfig.shrinklist < 0), "square", 'linewidth', 3, 'Color', 'b');
    plot(find(SEconfig.shrinklist > 0), res(SEconfig.shrinklist > 0), "square", 'linewidth', 3, 'Color', 'r');

    title(MatNames(fileNo));
    legend("ISI", "ISI with shrink", "shrink point", "enlarge point");
    xlabel("Iterations");
    ylabel("Res");
    set(gca,'FontSize',16);

    % save figure
    exportgraphics(gca, FigureOut(fileNo));
    hold off;

    % save data
    dlmwrite(DataOut, timeL(fileNo, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite(DataOut, iterL(fileNo, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite(DataOut, resL{fileNo}(1, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite(DataOut, resL{fileNo}(2, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite(DataOut, shrinklistL(fileNo, :), '-append', 'delimiter', ',', 'precision', 4);

end