%% Grit Specifier
% This module assigns the geometrical specification to each present grit 

function [gmin, gmax, GritSpec, time] = GSpecifier(gm, rmu, rsig, gmu, gsig, cmu, csig, wmu, wsig, hmu, hsig, nmu, nsig, omu, osig, smu, ssig)
    t = tic;
    % Counts the number of present grits
    k = 0;
    for i = 1 : 1 : size(gm, 1)
        for j = 1 : 1 : size(gm, 2)
            if gm(i, j) == 1
                k = k + 1;
            end
        end
    end
    
    % Creates a 8 x k matrix
    GritSpec = zeros(size(gm, 1), size(gm, 2), 8);
    
    n = sum(gm(gm == 1));
 
    rdist = makedist('Normal', 'mu', rmu, 'sigma', rsig);
    gdist = makedist('Normal', 'mu', gmu, 'sigma', gsig);
    cdist = makedist('Normal', 'mu', cmu, 'sigma', csig);
    wdist = makedist('Normal', 'mu', wmu, 'sigma', wsig);
    hdist = makedist('Normal', 'mu', hmu, 'sigma', hsig);
    ndist = makedist('Normal', 'mu', nmu, 'sigma', nsig);
    odist = makedist('Normal', 'mu', omu, 'sigma', osig);
    sdist = makedist('Normal', 'mu', smu, 'sigma', ssig);

    rrandom = rdist.random(n, 1);
    grandom = gdist.random(n, 1);
    crandom = cdist.random(n, 1);
    wrandom = wdist.random(n, 1);
    hrandom = hdist.random(n, 1);
    nrandom = ndist.random(n, 1);
    orandom = odist.random(n, 1);
    srandom = sdist.random(n, 1);
   
    pg = find(gm == 1);
    
    for i = 1 : 1 : size(gm, 1)
        for j = 1 : 1 : size(gm, 2)
            if gm(i, j) == 1
                idx1 = sub2ind(size(gm), i, j);
                idx = find(pg == idx1); 
                % row 1 = rake angle [degree]
                GritSpec(i, j, 1) = rrandom(idx);
                
                % row 2 = height of the grit [um]
                GritSpec(i, j, 2) = abs(grandom(idx));
                
                % row 3 = clearance angle [degree]
                GritSpec(i, j, 3) = abs(crandom(idx));
                
                % row 4 = width of the cut [mm]
                GritSpec(i, j, 4) = abs(wrandom(idx));
                
                % row 5 = hone radius [um]
                GritSpec(i, j, 5) = abs(hrandom(idx));
                
                % row 6 = nose radius [um]
                GritSpec(i, j, 6) = abs(nrandom(idx));

                % row 7 = oblique angle [degree]
                GritSpec(i, j, 7) = orandom(idx);
                
                % row 8 = stagnation angle [degree]
                GritSpec(i, j, 8) = srandom(idx);
            end
        end
    end

    GritSpec(GritSpec(:, :, 2) < 0, 2) = 0;
    GritSpec(GritSpec(:, :, 5) < 0, 5) = 0;
    gmin = max(min(GritSpec(:, :, 2), [], 'all'), 0);
    gmax = max(GritSpec(:, :, 2), [], 'all');
    time = ['The elapsed time for specifying the parameters is ', num2str(toc(t)), ' seconds'];

end
