 function [work, xfield, yfield, chipThickness, chipTotaleS, heattrack, p0, px, kt, kn, kr, f, fp, fc, fdmz, pplough, pcut, pdmz, ft, ftsmooth, fn, temp, tempmax, tempmax2, tempfield, epsi, q, time] = WSurface2(phii, xinter, zinter, xc0, zc0, xg, yg, zg, gmf, gm, gs, vc, tr, zerotemp, time, material, coolant, grit, ad, gc, pg, wd, wr, ww, wl, wh, werk, workdimz, divv, div, maxstep, rdoc, fdm, fdmethod, matcom, timeres, timeres2, Et, Rt, Ct, rhot, eps, SimDim, bcs, deltmax, plott, GST)
    t = tic;
    % Conversion from cm to mm
    ww = ww * 10;
    wd = wd * 10;
    pl = 0;
    syze1 = size(xg, 1);
    syze2 = size(xg, 2);
    syze3 = size(xg, 3);
    ag = 0;
    s = 1;
    tcounter = ag;
    % xdiv = linspace(0, max(xg(:, :, end), [], 'all') - (wr / 2.2), div);
    phii = rad2deg(phii);
    xdiv = linspace(0, workdimz(3), div);
    ydiv = linspace(0, workdimz(1), syze2);
    zdiv = linspace(0, workdimz(2), divv);
    [yfield, xfield, zfield] = meshgrid(ydiv, xdiv, zdiv);
    work = repmat(werk, [1, 1, maxstep]);
    tins = ones(size(xfield, 2), size(xfield, 1), size(xfield, 3)) .* zerotemp;
    tempfield = ones([size(tins), syze3]) .* zerotemp;
    qgen = ones(size(tins)) .* 1e-4;

    kt = zeros(size(xg));
    kn = zeros(size(xg));
    kr = zeros(size(xg));
    tau = zeros(size(xg));
    chipTotaleS = zeros(syze3, 1);
    
    temp = ones(syze3, 1) .* zerotemp;
    tempmax = ones(syze3, 1) .* max(zerotemp, [], 'all');
    q = zeros(syze3, 1);
    epsi = zeros(syze3, 1);
    Temp2Chip = 30;
    tchip = Temp2Chip;

    heatPattern = zeros(size(tins(:, :, 1)));
    heatWork = ones(size(tins)) .* .01;
    heattrack = zeros([size(tins(:, :, 1)), syze3]);
    AreaElements = 0;
    heatAreaPercent  = 0;
    support = heatAreaPercent;

    % fx = zeros(size(xg));
    % fy = zeros(size(xg));
    % fz = zeros(size(xg));

    f = zeros(syze3, 3);
    ft = ones(syze3, 1) * .1;
    ftsmooth = zeros(syze3, 1);
    fn = zeros(syze3, 1);
    fp = zeros(syze3, 3);
    fc = zeros(syze3, 3);
    fdmz = zeros(syze3, 3);

    pplough = zeros(syze3, 3);
    pcut = zeros(syze3, 3);
    pdmz = zeros(syze3, 3);
    
    % df = zeros(syze3, 3);
    chipThickness = cell(size(xg));
    chipThicknessRaw = cell(size(xg));
    InCutGrits = cell(syze3, syze2);
    cnt = [1, 1];
    theta = cell(size(xg));
    arcthet = zeros(size(xg));
    arclen = zeros(syze1, syze2);
    grds = wr + (gm .* 1e-3);
    qys = linspace(0, workdimz(3), div);
    Lc = sqrt(rdoc * wd);

    [Ny, Nx] = size(werk);
    [Nr, Nd, Nt] = size(xg);

    % Workpiece Initialist
    if any(xg(:, :, 1) >= -.01) & any(zg(:, :, 1) < wh) & s == 1
        xg_t = squeeze(xg(:, :, 1));
        zg_t = squeeze(zg(:, :, 1));

        indx = xg_t >= -.01;
        indz = zg_t <= wh + .01;
        [sqind(:, 1), sqind(:, 2)] = find(indx .* indz);
        diksz = unique(sqind(:, 2));

        for dd = 1 : numel(diksz)
            d = diksz(dd);
            rowrange = sqind(sqind(:, 2) == d);
            if numel(rowrange) == 1
                rowrange(end + 1) = rowrange + 1;
            end
            something(:, 1) = xg_t(min(rowrange) : max(rowrange), d);
            something(:, 2) = zg_t(min(rowrange) : max(rowrange), d);
            initialp1 = [something(1, 1), something(1, 2)];
            initialp2 = [something(end, 1), something(end, 2)];
                    % Grit trajectory index number increases CW
                    % Workpiece index number increases CCW
            [initialx, initialz, ~] = ArcPlotter(initialp1, initialp2, grds(min(rowrange), d), min(rowrange), d, s, ad, 'n');
            [~, puntoFinal] = min(abs(something(1, 1) - xdiv));
            [~, puntoDePartida] = min(abs(something(end, 1) - xdiv));
            workXComenzar = xdiv(puntoDePartida);
            workXFin = xdiv(puntoFinal);
            x_array = linspace(workXComenzar, workXFin, (puntoFinal - puntoDePartida + 1));
            z_interp = interp1(initialx, initialz, x_array, 'makima');
            z_interp(isnan(z_interp)) = wh;
            work(d, puntoDePartida : puntoFinal - 1, 1) = min(work(d, puntoDePartida : puntoFinal - 1, 1), z_interp(1 : end - 1));
            work(d, puntoDePartida : puntoFinal - 1, 2) = min(work(d, puntoDePartida : puntoFinal - 1, 2), z_interp(1 : end - 1));
            clear something
        end
    end
    % surf(work(:, :, 1))



    xz0 = [0, 0];
    dfl = zeros(1, 4);
    xdfl = 0;
    zdfl = 0;
    spec = struct( ...
      'E', Et, 'rho', rhot, 'A', pi * Rt ^ 2, 'L', wl, ...   % SI units
      'Ix', .25 * pi * Rt ^ 4, 'Iz', .25 * pi * Rt ^ 4, ...
      'zeta', [0.03 0.04], ...
      'k', (Et * pi * Rt ^ 2) / wl, ...
      'm', pi * wd * ww * rhot, ...
      'model_x','bending', 'G', 80e9, 'J', 2.4e-8, 'r_contact', 0.15, 'I_rot', 0.12, ...
      'm_lumped_x', 2.0, 'm_lumped_z', 6.0 );
    % initial state for vibration; everything set to zero other than k, c, and m
    state = struct('x', [0, 0], 'v', [0, 0], 'a', [0, 0]);

    for s = 2 : 1 : maxstep
        wcl = xc0 + fdm * time(s);
        [~, cnt(1)] = min(abs(wcl - qys));
        [~, cnt(2)] = min(abs((wcl + Lc) - qys));
        for d = 1 : 1 : syze2
            for r = pg{1, d}
                if xinter{r, d}(s) > 0 || xinter{r, d}(s - 1) > 0
                    [~, gind] = min(abs(xinter{r, d}(s) - xdiv));
                    if (zinter{r, d}(s - 1) < work(d, gind, s - 1) || zinter{r, d}(s) < work(d, gind, s - 1)) && (zinter{r, d}(s - 1) < wh || zinter{r, d}(s) < wh)
                        ag = ag + 1;
                        p1 = [xinter{r, d}(s - 1), zinter{r, d}(s - 1)];
                        p2 = [xinter{r, d}(s), zinter{r, d}(s)];
                        
                        [x, z, thet] = ArcPlotter(p1, p2, grds(r, d), r, d, s, ad, plott);

                        theta{r, d, s} = -thet;
                        arcthet(r, d, s) = abs(thet(1) - thet(end));
                        
                        pl = pl + 1;
                        chipThickness{r, d, s} = (work(d, gind) - z) .* abs(cosd(90 + thet));
                        [tsval, tsind] = min(abs(x(1) - xdiv));
                        [teval, teind] = min(abs(x(ad) - xdiv));
                        wrkstrtind = min(tsind, teind);
                        wrkendind = max(tsind, teind);
                        wz = interp1(x, z, xdiv(wrkstrtind : wrkendind), 'makima');
                        comparison = wz < work(d, wrkstrtind : wrkendind, s - 1);
                        work(d, wrkstrtind : wrkendind, s) = comparison .* wz + (1 - comparison) .* work(d, wrkstrtind : wrkendind, s - 1);
                        work(d, wrkstrtind : wrkendind, s + 1) = comparison .* wz + (1 - comparison) .* work(d, wrkstrtind : wrkendind, s - 1);
                        % surf(work(:, :, s))
                        % drawnow
                        chipThicknessRaw{r, d, s} = chipThickness{r, d, s};
                        chipThickness{r, d, s}(sign(chipThickness{r, d, s}) < 0) = 0;
                        InCutGrits{s, d}(end + 1) = r;
                        chipTotaleS(s, 1) = chipTotaleS(s, 1) + sum(chipThickness{r, d, s}, 'all');

                        [~, miin] = min(abs(min(-thet) - phii(:, s)));
                        [~, maax] = min(abs(max(-thet) - phii(:, s)));
                        miin = min(miin, maax);
                        maax = max(miin, maax);
                        heatPattern = gmf(miin : maax, :);
                        heatPattern = (imbinarize(imresize(heatPattern, [(abs(cnt(2) - cnt(1)) + 1), size(heatPattern, 2)]), "adaptive", "Sensitivity", 1))';
                        heatWork = ones(size(tins)) .* .01;
                        heatWork(:, cnt(1) : cnt(2), 1) = heatPattern;
                        heatWork(heatWork == 0) = .01;
                        heattrack(:, :, s) = heatWork(:, :, 1);
                        AreaElements = sum(heatWork(heatWork == 1), "all");
                        heatAreaPercent  = AreaElements / ((cnt(2) - cnt(1) + 1) * syze2);
                        support = heatAreaPercent;
                    end
                end
            end
        end
        work(:, :, s + 1) = work(:, :, s);


        if GST(4) == 1 && sum(~cellfun('isempty', InCutGrits(s, :))) > 0
            % GForce1 calculates the force coefficients
            arclen = grds .* deg2rad(arcthet(:, :, s));
            [kt1, kn1, kr1, phi_n, tau1, Vchip, mu, beta_n, p0, px, h_cr, signz, tau_f, l1, l2, l3, deltat{4, 1}] = GForce1(gm, gs, chipThickness(:, :, s), InCutGrits(s, :), vc, tr, temp(s), material);
            tau(:, :, s) = tau1;
            kt(:, :, s) = kt1;
            kn(:, :, s) = kn1;
            kr(:, :, s) = kr1;
            % GForce2 calculaates forces for each time-step
            [f(s, 1), f(s, 2), f(s, 3), fp(s, 1), fc(s, 1), fdmz(s, 1), fp(s, 3), fc(s, 3), fdmz(s, 3), pplough(s, :), pcut(s, :), pdmz(s, :), ft(s, 1), fn(s, 1), thm, deltat{5, 1}] = GForce2(ad, h_cr, chipThickness(:, :, s), InCutGrits(s, :), p0, px, tau(:, :, s), signz, l1, l2, l3, mu, kt(:, :, s), kn(:, :, s), kr(:, :, s), gs, rdoc, theta(:, :, s), s);

            fti = ft(s, 1);
            if GST(5) == 1 && s > 360
                tfluid = tr - deltmax;
                [allin] = GTempFD(time(1, s), time(1, s + 1), timeres, wcl, workdimz, work(:, :, s), tins, heatWork, heatAreaPercent, div, xg(:, :, s), f(s, 3), tr, tfluid, material, fdmethod, coolant, grit, rdoc, wd, ww, wr, gm, vc, eps, bcs, deltmax, tchip, fdm);
                q(s) = mean(allin.BC.z1.neumann_params.q, 'all');
                [tempfield(:, :, :, s), deltat{6, 1}] = GTempFD3D3(allin, s);
                tempfield(:, :, :, s : end) = repmat(tempfield(:, :, :, s), 1, 1, 1, (maxstep - s + 2));
            end
            tins = tempfield(:, :, :, s);
            temp(s) = mean(tins(2 : end - 1, cnt, end), "all");
            tins2 = mean(tins(2 : end, cnt, end - 1), 'all');
            tempmax(s) = max(tins(2 : end - 1, cnt, end), [], "all");
            tempmax2(s) = tins2;
            tchip = max(tempmax(s), Temp2Chip);
        elseif GST(4) == 1 && sum(~cellfun('isempty', InCutGrits(s, :))) == 0 && GST(5) == 1 && s > 360
            tfluid = tr - deltmax;
            [allin] = GTempFD(time(1, s), time(1, s + 1), timeres, wcl, workdimz, work(:, :, s), tins, heatWork, heatAreaPercent, div, xg(:, :, s), f(s, 3), tr, tfluid, material, fdmethod, coolant, grit, rdoc, wd, ww, wr, gm, vc, eps, bcs, deltmax, tchip, fdm);
            q(s) = mean(allin.BC.z1.neumann_params.q, 'all');
            [tempfield(:, :, :, s), deltat{6, 1}] = GTempFD3D3(allin, s);
            tempfield(:, :, :, s : end) = repmat(tempfield(:, :, :, s), 1, 1, 1, (maxstep - s + 2));
            tins = tempfield(:, :, :, s);
            temp(s) = mean(tins(2 : end - 1, cnt, end), "all");
            tins2 = mean(tins(2 : end, cnt, end - 1), 'all');
            tempmax(s) = max(tins(2 : end - 1, cnt, end), [], "all");
            tempmax2(s) = tins2;
            tchip = max(tempmax(s), Temp2Chip);
            % temp = mean(tins(:, cnt (1) : cnt(2)), "all");
        end
        if GST(6) == 1
            [vibparams] = preDeflector(spec);
            [state, history, deltat{7, 1}] = Deflector(state, [ft(s, 1), fn(s, 1)], support, (time(1, s + 1) - time(1, s)), vibparams);
            % xdfl and zdfl are the relative deflection of the tool
            % between two consequtive time steps,
            % xdfl = xz(1) - xz0(1);
            % zdfl = xz(3) - xz0(2);
            % xinter = cellfun(@(x) x + xdfl, xinter, 'UniformOutput', false);
            % zinter = cellfun(@(z) z - zdfl, zinter, 'UniformOutput', false);
            % xg1 = xg + xdfl;
            % zg1 = zg + zdfl;
            % dfl = xz;
            % xz0 = xz;
        end

        %figure(1)
        %surf(work(:, :, s))    
        %surf(work(:, 1:100, s))
        %view(-360, 0)
        %drawnow
        %tins = tempfield;
    end
    time = ['The elapsed time for computing the worked surface is ', num2str(toc(t)), ' seconds'];
end