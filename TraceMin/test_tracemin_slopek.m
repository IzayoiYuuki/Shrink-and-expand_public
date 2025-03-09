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

OutNames = './Figure/' + MatNames + '_TM_slopek.pdf';

Nevs = [100; 100; 100; 100; 104; 156; 199; 384; 481; 500; 500; 500];
Maxiter = 500;
tol = 1e-10;

Sconfig = struct();
Sconfig.method = 'pcg';
Sconfig.res = 1e-5;
Sconfig.iter = 5;

SEconfig = struct();
SEconfig.rule = 'slopek';
SEconfig.slopestep = 10;
SEconfig.enlargetol = 1.1;
SEconfig.enlargesteps = 2;
warmupiter = 5;
SEconfig.warmuptol = 1e-4;

dlmwrite('./Figure/Data_TM_slopek.txt', date, '-append', 'delimiter', '', 'precision', 4);

for fileNo = 1 : 8

    % load matrix
    disp(MatNames(fileNo));
    dlmwrite('./Figure/Data_TM_slopek.txt', fileNo, '-append', 'delimiter', '', 'precision', 4);

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
    
    % Without shrink
    SEconfig.warmupiter = Maxiter;
    ftime = tic;
    [~, ~, iter, res] = myTraceMin(A, B, X, nev, tol, Maxiter, Sconfig, SEconfig);
    timeL(fileNo, 1) = toc(ftime);
    iterL(fileNo, 1) = iter;
    resL{fileNo}(1, 1:length(res)) = res;
    
    semilogy((1:iter), res, '-o', 'linewidth', 2);
    hold on;
    
    % With shrink
    SEconfig.warmupiter = warmupiter;
    ftime = tic;
    [~, ~, iter, res, SEconfig.shrinklist] = myTraceMin(A, B, X, nev, tol, Maxiter, Sconfig, SEconfig);
    timeL(fileNo, 2) = toc(ftime);
    iterL(fileNo, 2) = iter;
    resL{fileNo}(2, 1:length(res)) = res;
    shrinklistL(fileNo, 1:length(SEconfig.shrinklist)) = SEconfig.shrinklist;
    
    semilogy((1:iter), res, '-o', 'linewidth', 2); hold on;
    SEconfig.shrinklist = SEconfig.shrinklist(1:iter);
    plot(find(SEconfig.shrinklist < 0), res(SEconfig.shrinklist < 0), "square", 'linewidth', 3, 'Color', 'b');
    plot(find(SEconfig.shrinklist > 0), res(SEconfig.shrinklist > 0), "square", 'linewidth', 3, 'Color', 'r');
    
    title(MatNames(fileNo));
    legend("TraceMin", "TraceMin with shrink", "shrink point", "enlarge point");
    xlabel("Iterations");
    ylabel("Res");
    set(gca,'FontSize',16);
    
    % save figure
    exportgraphics(gca, OutNames(fileNo));
    hold off;
    
    % save data
    dlmwrite('./Figure/Data_TM_slopek.txt', timeL(fileNo, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Figure/Data_TM_slopek.txt', iterL(fileNo, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Figure/Data_TM_slopek.txt', resL{fileNo}(1, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Figure/Data_TM_slopek.txt', resL{fileNo}(2, :), '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite('./Figure/Data_TM_slopek.txt', shrinklistL(fileNo, :), '-append', 'delimiter', ',', 'precision', 4);
    
end