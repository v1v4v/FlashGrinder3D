
function [k, presentIndex, pa, phis, phie, ftr, work, xg0, zg0, time] = GParameters(gm, xc0, zc0, rdoc, rds, ww, fdr, wh, gmu, remRatio, om, div, ro)
    t = tic;
    % k is the number of cutting edges in each disk
    k = zeros(1, size(gm, 2));
    for j = 1 : 1 : size(gm, 2)
        presentIndex{j} = find(gm(:, j) ~= 0)';
        %presentIndex{j}(end + 1) = presentIndex{j}(1);
        for i = 1 : 1 : size(gm, 1)
            if gm(i,j) ~= 0
                k(j) = k(j) + 1;
            end
        end
    end
    ww = ww * 10;
    % Static active grits
    StaticActiveGrits = sum(k, 'all') / (ww * 2 * rds * pi);
    StaticActiveGrits2 = sum(k, 'all') / (ww * 2 * rds * pi * remRatio);
    disp(['Static number of active grits is ', num2str(StaticActiveGrits), ' grits / mm^2'])
    disp(['Static number of active grits (before patterning) is ', num2str(StaticActiveGrits2), ' grits / mm^2'])    

    % Pitch angle [degree]
    pa = 360 ./ k;
    % Feed per grit(tooth) per revolution [mm / rev]
    ftr = fdr ./ k;

    % Start and exit angle
    for c = 1 : 1 : size(gm, 2)
        Rmax = rds + max(gm(:, c), [], 'all') * 1e-3;
        for p = 2 : 1 : length(presentIndex{1, c})
            n = presentIndex{1, c}(p);
            m = presentIndex{1, c}(p - 1);
            if om == "down"
                phis(n, c) = pi - acos(1 - ((rdoc - Rmax + (rds + (gm(n, c) * 1e-3))) / (rds + (gm(n, c) * 1e-3))));
                phie(n, c) = pi - asin(((gm(n, c) - gm(m, c)) * 1e-3) / ((n - m) * ftr(c)));
            elseif om == "up"
                phis(n, c) = asin(((gm(n, c) - gm(m, c)) * 1e-3) / ((n - m) * ftr(c)));
                phie(n, c) = acos(1 - ((rdoc - Rmax + (rds + (gm(n, c) * 1e-3))) / (rds + (gm(n, c) * 1e-3))));
            end
        end
    end
    %%
    workdist = makedist("HalfNormal", "mu", wh, "sigma", 0);

    work = workdist.random(size(gm, 2), div);

    work = 2 * wh - work;
    
    % Initial position of the peripheral elements

    % fi = (linspace(0, -2 * pi, size(gm, 1)))';
    % fi = repmat(fi, 1, size(gm, 2));
    % 
    % xg0 = xc0 + (((gm .* 1e-3) + rds) .* cos(fi));
    % zg0 = (zc0 + rds + wh - rdoc + (gmax * 1e-3)) + (((gm .* 1e-3) + rds) .* sin(fi));

    % xg0 = xc0 .* ones(size(gm));
    % zg0 = rds + wh + zc0 + (gmax .* 1e-3) - rdoc;
    % gmu = max(gm, [], "all");
    xg0 = zeros(size(gm, 1), 0);
    zg0 = zeros(size(gm, 1), 0);
    for f = 1 : 1 : size(gm, 1)
        fi = - f * ((2 * pi) / size(gm, 1)) + (3 * pi / 2);
        xg0(f) = xc0;
        zg0(f) = rds + wh + zc0 + ro(2) + (gmu * 1e-3) - rdoc;
    end

    time = ['The elapsed time for parameter function is ', num2str(toc(t)), ' seconds'];
end