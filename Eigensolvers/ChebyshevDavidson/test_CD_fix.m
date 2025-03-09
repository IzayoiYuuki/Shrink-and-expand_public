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
FigureOut = fullfile(RootPath, "Figure", MatNames + "_CD_fix.pdf");
DataOut = fullfile(RootPath, "Figure", "Data_CD_fix.txt");

% Problem setting
Nevs = [100; 100; 100; 100; 104; 156; 199; 384; 481; 500; 500; 500];
Maxiter = 1000;
tol = 1e-10;

% solver setting
SEconfig.rule = 'fix';
SEconfig.enlargesteps = 2;
SEconfig.shrinksteps = 10;
warmupiter = 30;
SEconfig.warmuptol = 1;

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
    nex = ceil(nev/4);

    % initial guess
    rng(0);
    [n, ~] = size(A);
    X = randn(n, nex);

    % parameters of Chebyshev-Davidson
    et = eigs(A, nev + 10, 'smallestabs');
    CDconfig.lowb = et(nev + 10);
    CDconfig.upb = norm(A, 1);
    CDconfig.polyorder = 25;
    CDconfig.submax = ceil(1.5*nev);
    CDconfig.newsub = nex;

    % CD without shrink
    SEconfig.warmupiter = Maxiter;
    tic;
    [~, ~, iter, res, logs] =...
        myChebyshevDavidson(A, X, nev, tol, Maxiter, CDconfig, SEconfig);
    log_restart = logs.restart;
    res_now_1 = logs.res_now;
    nconv_1 = logs.nconv;
    timeL(fileNo, 1) = toc;
    iterL(fileNo, 1) = iter;
    resL{fileNo}(1, 1:length(res)) = res;

    semilogy((1:iter), res, '-*', 'linewidth', 2);
    hold on;

    % CD with shrink
    SEconfig.warmupiter = warmupiter;
    tic;
    [~, ~, iter, res, logs] =...
        myChebyshevDavidson(A, X, nev, tol, Maxiter, CDconfig, SEconfig);
    log_restart = logs.restart;
    res_now_2 = logs.res_now;
    nconv_2 = logs.nconv;
    SEconfig.shrinklist = logs.shrinklist;
    timeL(fileNo, 2) = toc;
    iterL(fileNo, 2) = iter;
    resL{fileNo}(2, 1:length(res)) = res;
    shrinklistL(fileNo, 1:length(SEconfig.shrinklist)) = SEconfig.shrinklist;

    semilogy((1:iter), res, '-o', 'linewidth', 2); hold on;
    SEconfig.shrinklist = SEconfig.shrinklist(1:iter);
    plot(find(SEconfig.shrinklist < 0), res(SEconfig.shrinklist < 0), "square", 'linewidth', 3, 'Color', 'b');
    plot(find(SEconfig.shrinklist > 0), res(SEconfig.shrinklist > 0), "square", 'linewidth', 3, 'Color', 'r');

    title(MatNames(fileNo));
    legend("Chebyshev-Davidson", "Chebyshev-Davidson with shrink", "shrink point", "enlarge point");
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