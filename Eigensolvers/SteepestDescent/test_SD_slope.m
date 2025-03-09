clear;
hold off;

WarmUp();

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

FileNames = './Matrices/' + MatNames + '.mat';

OutNames = './Figure/' + MatNames + '_SD_slope.pdf';

Nevs = [100; 100; 100; 100; 104; 156; 199; 384; 481; 500; 500; 500];
Maxiter = 3500;
tol = 1e-10;

SEconfig.rule = 'slope';
SEconfig.enlargetol = 1.1;
SEconfig.enlargesteps = 2;
warmupiter = 5;
SEconfig.warmuptol = 1e-4;

dlmwrite('./Figure/Data_SD_slope.txt', date, '-append', 'delimiter', '', 'precision', 4);
for fileNo = 1 : 11

    disp(MatNames(fileNo));
    dlmwrite('./Figure/Data_SD_slope.txt', fileNo, '-append', 'delimiter', '', 'precision', 4);

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

    % LOBPCG without shrink
    SEconfig.warmupiter = Maxiter;
    tic;
    [~, ~, iter, res, ~] = mySteepestDescent(A, X, Myprec, nev, tol, Maxiter, SEconfig);
    timeL(fileNo, 1) = toc;
    iterL(fileNo, 1) = iter;
    resL{fileNo}(1, 1:length(res)) = res;

    semilogy((1:iter), res, '-*', 'linewidth', 2);
    hold on;

    % LOBPCG with shrink
    SEconfig.warmupiter = warmupiter;
    tic;
    [~, ~, iter, res, SEconfig.shrinklist] = mySteepestDescent(A, X, Myprec, nev, tol, Maxiter, SEconfig);
    timeL(fileNo, 2) = toc;
    iterL(fileNo, 2) = iter;
    resL{fileNo}(2, 1:length(res)) = res;
    shrinklistL(fileNo, 1:length(SEconfig.shrinklist)) = SEconfig.shrinklist;

    semilogy((1:iter), res, '-o', 'linewidth', 2); hold on;
    SEconfig.shrinklist = SEconfig.shrinklist(1:iter);
    plot(find(SEconfig.shrinklist < 0), res(SEconfig.shrinklist < 0), "square", 'linewidth', 3, 'Color', 'b');
    plot(find(SEconfig.shrinklist > 0), res(SEconfig.shrinklist > 0), "square", 'linewidth', 3, 'Color', 'r');

    hold on;

    legend("SD", "SD with shrink", "shrink point", "enlarge point");
    xlabel("Iterations");
    ylabel("Res");
    set(gca,'FontSize',16);

    % save figure
    exportgraphics(gca, OutNames(fileNo));
    hold off;

    % save data
    dlmwrite('./Figure/Data_SD_slope.txt', timeL(fileNo, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Figure/Data_SD_slope.txt', iterL(fileNo, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Figure/Data_SD_slope.txt', resL{fileNo}(1, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Figure/Data_SD_slope.txt', resL{fileNo}(2, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Figure/Data_SD_slope.txt', shrinklistL(fileNo, :), '-append', 'delimiter', ',', 'precision', 4);

end
