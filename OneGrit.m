
warning('off')
maxNumCompThreads('automatic');
disp('Grinding Gedanken Twin')
disp('Developed at Manufacturing Research Laboratory (MRL), Sabanci University, 2.022e3')
disp('All rights reserved \copyright')
disp('                    ')

disp('Boundary condition numbering for temperature solver')

disp('            x1              ')
disp('   ___O~~~~~~~~~~~~\        ')
disp('  | \___O~~~~~~~~~~~\       ')
disp('  \  \___O~~~~~~z1~~~\      ')
disp('   \  \___O~~~WP~~~~~~\   y1')
disp(' y0 \  \___O~~~~~~~~~~~\    ')
disp('     \  \___O~~~~~~~~~~~\   ')
disp('      \ |    |   x0     |   ')
disp('       \|____|__________|   ')
disp('                    z0      ')
disp('      >>----feed dir.---->> ')

%% Input

% General simulation tweaks (GST)
% 1 -> Patterned conventional wheel (1) or not(0)
% 2 -> Point Grit (1) or shaped grit (0)
% 3 -> Create movie (1) or not (0) / NOTE: Higher precision is available when movie is off.
% 4 -> Force (1)
% 5 -> Temperature (1)
% 6 -> Beam Deflection on (1)
% 7 -> Wear on (1)

% GST = [0, 1, 2, 0, 0];
%
% Specific tweaks
%
% WD = 2; %  Wheel Diameter [cm]
% WR = (WD * 10) / 2; %  Wheel radius [mm]
% WW = 1; %  Wheel Width [cm]
% WH = 25; %  Work Height [mm]
% SN = 5; %  Number of strips in the width of the wheel [#]
% SA = 79; %  Strip Angle [degrees]
% Fn1 = 10; %  Fineness or Number of divisions between two adjacent strips
% Fn2 = 20; %  Fineness in periphery of the wheel
% NR = 5; %  Number of rotations
% TS = 'Dgr'; %  Trigonometry system, "Rad" or "Dgr"
% MaxStep = 1200; %  The maximum time step, must be less than NR * 360 / stpA
% EA = 360; %  End Angle or Last angle
% disks = SN * (Fn1 + 1); %  Number of elemental disks
% Pttrn = "/"; %  Pattern type, /, \, v, ^, x, or arbitrary curve ~ (the equation should be defined in 'EQ' variable)
% EQ = "(x ^ 5 + x ^ 3) / (x ^ 2 + 1)"; %  Equation of the pattern defined in x
% INT = [-1 1]; %  The interval for the curve
% PD = "hor"; %  Pattern direction, either 'ver'  for vertical or 'hor' for horizontal
% PTh = 7; %  Pattern thickness, 0, 1, 2, 3, 4, ...
% PS = 0; % Pattern spacing, white space between tow pattern rows
%
% GD = .1;  % Grit density, expressed as the percentage of the grits that initially are present
% GSh = "m"; %  Grit shape (^, /, \, -, ~, or m), if curved, the equation should be provided in the GEQ in x and y, or the grit shape can be input manually
% GR = [4, 3]; %  Grit resolution, odd value for the number of columns please
% GritTemplate = [1 3 4 5 5.5 5 4 3 1; .8 2.5 3.2 4 4.5 4 3.2 2.5 .8; .2 1.5 2 2 2.2 2 2 1.5 .2; 0 0 0 .1 .2 .1 0 0 0]; %  Grit shape, manually defined in a matrix
% GEQ = "cos(y) + sin(x)";
% xint = [-2 2]; %  x interval for the grit equation
% yint = [-4 4]; %  y interval for the grit equation
%
% GM = []; %  Grit Matrix
% GMin = 130; %  Minimum grit height in um
% GMax = 150; %  Maximum grit height in um
% rpm = 2000; %  Revolution per minute [rpm]
% rps = rpm / 60; %  Revolution per second [rps]
% Omega = rps * 2 * pi; %  Rotational speed [rad/s]
% StpA = 1; %  Step angle or Angle increment, current state in Caps., DEGREE or radian
% fdm = 50; %  Feed [mm/s]
% %fdmm = Omega * WR; % Rolling speed of the wheel, equal to the tangential speed of the rolling (w/o sliding) wheel [mm/s]
% %fdm = 20000 + (-fdmm * 60); %  Feed per minute, (-fdmm * 60) is the inherent rolling speed which is subtracted [mm/min]
% fdr = fdm / rps; %  Feed per revolution [mm/rev]
% OM = 'down'; %  Operation mode, down/up grinding
% RDOC = .1; %  Radial depth of cut [mm]
% ADOC = WW * 10; %  Axial depth of cut, or the width of the wheel [mm]
% RO = .01; %  Run Out [um]
% phi = 30; %  Rotation of the wheel

if strcmp(wheel, 'create')
    %% Pattern generation
    % GPatterner function create the pattern based on the input, GM is the full wheel surface
    [GM, GMF, remRatio, deltat{1, 1}] = GPatterner2(WD, WW, SN, SA, Fn1, Fn2, Pttrn, EQ, INT, PD, PTh, symm, PS, GD, wheel, GST);
    
    %% Grit specifications
    % GSpecifier function assigns the certain geometrical parameters to each present grit
    [gmin, gmax, GS, deltat{2, 1}] = GSpecifier(GM, rmu, rsig, gmu, gsig, cmu, csig, wmu, wsig, hmu, hsig, nmu, nsig, omu, osig, smu, ssig);
elseif strcmp(wheel, 'load')
    load([MyDir, '\GS.mat'])
    GM = GS(:, :, 2);
    GMF = GS(:, :, 2);
    gmin = max(min(GS(:, :, 2), [], 'all'), 0);
    gmax = max(GS(:, :, 2), [], 'all');
elseif strcmp(wheel, 'loadNpattern')
    load([MyDir, '\GS.mat'])
    [GM, GMF, remRatio, deltat{1, 1}] = GPatterner2(WD, WW, SN, SA, Fn1, Fn2, Pttrn, EQ, INT, PD, PTh, symm, PS, GD, wheel, GST);
    gmin = max(min(GS(:, :, 2), [], 'all'), 0);
    gmax = max(GS(:, :, 2), [], 'all');
    GS2 = zeros([size(GM), size(GS, 3)]);
    for gss = 1 : size(GS, 3) 
        GS2(:, :, gss) = GS(1 : size(GM, 1), 1 : size(GM, 2), gss);
    end
    GS = GS2 .* GM;
    clear GS2
end

%% Grit shaper
% GShaper creates the grit matrix(GM) using specified shaped grits
[GM, GritTemp, GSS, deltat{3, 1}] = GShaper(GM, GS(:, :, 2), GSh, GR, xint, yint, GritTemplate, GEQ, gmin, gmax, GST);

% GM = distributed(GM);
% cutting velocity in [mm/s]
vc = 2 * pi * (WR + (GM * 1e-3)) * rps;

%% Grit Force
% This fuction takes in the chip thicknesses and spits the forces in three
% directions out
% if GST(4)
%     [tau1, kt, kn, kr, phi_n, tau, Vchip, mu, beta_n, p0, h_cr, tau_f, deltat{4, 1}] = GForce1(GM, GS, vc, Tr, material);
% end

%% Some parameter calculation
[k, PG, PA, phis, phie, ftr, Work, xg0, zg0, deltat{5, 1}] = GParameters(GM, xc0, zc0, RDOC, WR, WW, fdr, WH, gmu, remRatio, OM, div2, RO);

%% Grit Tracker
% GTracker function predicts the trajectory of each grit and element
[xg, yg, zg, phi, phii, time, deltat{6, 1}] = GTracker4(GM, WR, fdm, Omega, StpA, NR, xg0, zg0, VirChat, exen, RO, phase, offcentre, OM);
MaxStep
% [xg2, zg2, ~, ~, ~] = GTracker4(GH, WR, fdm, Omega, StpA, NR, RO, xg0, zg0, OM);

%% Trajectory Intersections
% "trim" variable inside the function must be set accordingly, otherwise the code will crash.
% it must be greatest value less than approximately 1/4 of the number of rows in the gm matrix
% fit model can be found from matlab help, 'linearinterp' works better if possible, if not 'poly9'
[xinter, zinter, deltat{7, 1}] = WIntersections(xg, zg, GM, PG, WH, gmax, fitModel);

%% Main engine, Surface roughness, Force, and Temperature
[worked, workXfield, workYfield, chth, chipTotale, heattrack, p0, px, kt, kn, kr, f, fp, fc, fdmz, pplough, pcut, pdmz, ft, ftsmooth, fn, temp, tempmax, tempmax2, tempfield, epsiloan, q, deltat{8, 1}] = WSurface2(phii, xinter, zinter, xc0, zc0, xg, yg, zg, GMF, GM, GS, vc, Tr, Tins, time, material, coolant, grit, ad, GC, PG, WD, WR, WW, WL, WH, Work, workDimz, div1, div2, MaxStep, RDOC, fdm, FDmethod, MatCom, timeres, timeres2, Et, Rt, Ct, rhot, eps, SimDim, BC, deltmax, plott, GST);
f = -f;
wxf1 = workXfield(:, :, 1)';
wyf1 = workYfield(:, :, 1)';
% Removed material volume [mm^3]
RemovedMaterial = sum(WH - worked(:, :, end), 'all') * (workDimz(1) / size(GM, 2)) * (workDimz(3) / div2);
% MRR [mm^3 / s]
MRR = RemovedMaterial / time(end);

% Average chip thickness
chipMax = cellfun(@(chip) max(chip), chth, 'UniformOutput', false);
chipEmpty = cellfun(@isempty, chipMax);
chipMaxMean = nan(size(chipEmpty));
chipMaxMean(~chipEmpty) = cellfun(@(chipp) chipp, chipMax(~chipEmpty));
clear chipMax
chipMean = mean(chipMaxMean, 3, "omitnan");
clear chipMaxMean
chipMean(isnan(chipMean)) = 0;
chipMeanNZ = chipMean(chipMean > 0);
chipMeanVector = reshape(chipMeanNZ, [numel(chipMeanNZ) 1]);
chipdist = fitdist(chipMeanVector, 'Rayleigh');


% Surface properties
% This function calculates the surface roughness in different approaches
[RAx, RAy, RMSx, RMSy, deltat{9, 1}] = SRoughness(worked(:, :, end), xinter, workXfield(:, 1, 1));
RaAve = [mean(RAx, "all"), mean(RAy, "all")];
disp("RA, x 'n y")
disp(RaAve(1))
disp(RaAve(2))
disp('')
RMSAve = [mean(RMSx, "all"), mean(RMSy, "all")];
disp("RMS, x 'n y")
disp(RMSAve(1))
disp(RMSAve(2))

fxx = f(361 : end, 1);
fzz = f(361 : end, 3);
fmean(1) = mean(fxx, 'all');
fmean(3) = mean(fzz, 'all');
fmeanuq(1) = mean(fxx(fxx > fmean(1)), 'all');
fmeanuq(3) = mean(fzz(fzz < fmean(3)), 'all');
fmax(1) = min(fxx);
fmax(3) = max(fzz);
fprintf('<strong> Mean force [N/mm] </strong> \n')
disp(fmean ./ workDimz(1))

fxp = fp(361 : end, 1);
fzp = fp(361 : end, 3);
fpmean(1) = mean(fxp, 'all');
fpmean(3) = mean(fzp, 'all');
fprintf('<strong> Mean ploughing force [N/mm] </strong> \n')
disp(fpmean ./ workDimz(1))
disp('Ploughing force percentage')
% disp(mean(pplough, "omitnan") .* 100)
disp((fpmean ./ fmean) .* 100)

fxdmz = fdmz(361 : end, 1);
fzdmz = fdmz(361 : end, 3);
fdmzmean(1) = mean(fxdmz, 'all');
fdmzmean(3) = mean(fzdmz, 'all');
fprintf('<strong> Mean DMZ force [N/mm] </strong> \n')
disp(fdmzmean ./ workDimz(1))
disp('DMZ force percentage')
% disp(mean(pdmz, "omitnan") .* 100)
disp((fdmzmean ./ fmean) .* 100)

fxc = fc(361 : end, 1);
fzc = fc(361 : end, 3);
fcmean(1) = mean(fxc, 'all');
fcmean(3) = mean(fzc, 'all');
fprintf('<strong> Mean cutting force [N/mm] </strong> \n')
disp(fcmean ./ workDimz(1))
disp('Cutting force percentage')
% disp(mean(pcut, "omitnan") .* 100)
disp((fcmean ./ fmean) .* 100)
if strcmp(chipdist.DistributionName, "Normal")
    disp('Overalll undeformed chip average and STD: ')
    disp([chipdist.mu chipdist.sigma])
else
    disp('Rayleigh distribution scale parameter:')
    disp([chipdist.B])
end

save([MyDir, '\testResult', num2str(TestNo), '.mat'], 'GS', 'worked', 'RaAve', 'RMSAve', 'RemovedMaterial', 'MRR', 'chipdist', 'chipMean', 'fpmean', 'fdmzmean', 'fcmean', 'pplough', 'pcut', 'pdmz', 'ft', 'ftsmooth', 'fn', 'f', 'fmean', 'fmeanuq', 'fmax', 'tempfield', 'temp', 'tempmax')


%% Temperature
%[Tx, Ty, Tz, deltat{10, 1}] = GTemp(~, ~);

%% Exhilarater
% This section is for making the movie
% figure(1)
% set(figure(1), 'Position', get(0, 'Screensize'));
% movieVector(s) = getframe(figure(1));
% movieWriter = VideoWriter('chip', 'MPEG-4');
% myWriter.FrameRate = 30;
% open(movieWriter);
% writeVideo(movieWriter, movieVector2);
% close(movieWriter);

%for j = 1 : 720; for i = 1 : 90; plot3(xg(:, i, j), zg(:, i, j), ycartz(i, :)); hold on; grid minor; end; drawnow; clf;end
%for s = 1 : 3600;figure(1); grid minor;plot3(xg(:, :, s), yg, zg(:, :, s)); drawnow; end

beep on