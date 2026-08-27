%%  Grit Tracker 4
%   Instantaneous GTracker
%   'specs' is just the grit height vector

function [xg, yg, zg, phi, phii, tt, time] = GTracker4(gm, rds, ft, omg, stpa, nr, xg0, zg0, virchat, exen, ro, phase, offcentre, om)
    t = tic;
    % nr * 360 is the total number of steps, nr rotation(s) times the number of steps
    % in each rotation which is 360 (each step = 1 deg)

    % offcentre and ro are deemed the same
    
    gm1 = size(gm, 1);
    gm2 = size(gm, 2);
    MaxStepTraj = round(nr * 360)
    MaxTraj = MaxStepTraj / stpa
    phii = linspace(0, 2 * pi, gm1);
    phii = phii';
    ppp = linspace(0, nr * 2 * pi, round(nr * 360 * (1 / stpa)));
    tt = (ppp ./ omg);
    ppp = phii + ppp;
    phi = zeros(gm1, gm2, round(nr * 360 * (1 / stpa)));
    xg = zeros(gm1, gm2, round(nr * 360 * (1 / stpa)));
    zg = zeros(gm1, gm2, round(nr * 360 * (1 / stpa)));
    
    
    % Conversion to mm/s
    %ft = ft / 60;

    phii = repmat(phii, 1, MaxStepTraj);
    phii = cell2mat(arrayfun(@(c) circshift(phii(:,c), c-1), 1:size(phii,2), 'UniformOutput', false));
    
    for qq = 1 : 1 : size(phi, 3)
        phi(:, :, qq) = repmat(ppp(:, qq), 1, gm2);
    end

    phi2 = rem(phi, 2 * pi);
    
    % Calculation the time of each step

    xg0 = xg0';
    zg0 = zg0';

    xg0 = repmat(xg0, 1, gm2);    
    zg0 = repmat(zg0, 1, gm2);

    % Virtual Chatter
    
    
    % This part takes care of eccentricity of the wheel (if exists)
    thet = linspace(0, 2 * pi, gm1);
    mx = max((1/exen) * rds, exen * rds);
    mn = min((1/exen) * rds, exen * rds);
    exenmat = (mx * mn) ./ sqrt((mn .* cos(thet)) .^ 2 + (mx .* sin(thet)) .^ 2);
    %exenmat = exenmat + max(exenmat);
    exenmat = circshift(exenmat', round((phase / (2 * pi)) * gm1));
    exenmat = repmat(exenmat, 1, gm2);

    % Run-Out
    freq = (2 * pi) / ro(3);
    ROVector = ro(1) * 1e3 + (ro(2) - ro(1)) * 1e3 * (sin(freq * thet) + 1) / 2;
    ROTemplate = repmat(ROVector', 1, gm2);
    ROTemplate(gm == 0) = 0;

    gm = gm + ROTemplate;
    
    % Loops over the steps
    % Operation mode, down or up milling
    if size(gm, 3) == 1
        xg00 = (exenmat(:, :) .* phi(:, :, 1) - (gm(:, :) .* 1e-3) .* sin(phi(:, : , 1)) + ft * tt(1));
        if om == "down"
            for s = 1 : 1 : round(nr * 360 * (1 / stpa))
                xg(:, :, s) = (exenmat(:, :) + gm(:, :) .* 1e-3) .* cos(phi(:, :, s)) + (ft * tt(s)) + xg0 + (offcentre * cos((-(s * stpa) * pi / 180) + phase));
                zg(:, :, s) = -(exenmat(:, :) + gm(:, :) .* 1e-3) .* sin(phi(:, :, s)) + zg0 + (offcentre * sin((-(s * stpa) * pi / 180) + phase));
            end
        elseif om == "up"
            for s = 1 : 1 : round(nr * 360 * (1 / stpa))
                xg(:, :, s) = (exenmat(:, :) + gm(:, :) .* 1e-3) .* cos(phi(:, :, s)) + ft * tt(s) + xg0 + (offcentre * cos((-(s * stpa) * pi / 180) + phase));
                zg(:, :, s) = (exenmat(:, :) + gm(:, :) .* 1e-3) .* sin(phi(:, :, s)) + zg0 + (offcentre * sin((-(s * stpa) * pi / 180) + phase));
            end
        end
    elseif size(gm, 3) > 1
        xg00 = (exenmat(:, :) .* phi(:, :, 1) - (gm(:, :, 1) .* 1e-3) .* sin(phi(:, :, 1)) + ft * tt(1));
        if om == "down"
            for s = 1 : 1 : round(nr * 360 * (1 / stpa))
                xg(:, :, s) = (exenmat(:, :) + gm(:, :, s) .* 1e-3) .* cos(phi(:, :, s)) + ft * tt(s) + xg0 + (offcentre * cos((-(s * stpa) * pi / 180) + phase));
                zg(:, :, s) = -(exenmat(:, :) + gm(:, :, s) .* 1e-3) .* sin(phi(:, :, s)) + zg0 + (offcentre * sin((-(s * stpa) * pi / 180) + phase));
            end
        elseif om == "up"
            for s = 1 : 1 : round(nr * 360 * (1 / stpa))
                xg(:, :, s) = (exenmat(:, :) + gm(:, :, s) .* 1e-3) .* cos(phi(:, :, s)) + ft * tt(s) + xg0 + (offcentre * cos((-(s * stpa) * pi / 180) + phase));
                zg(:, :, s) = (exenmat(:, :) + gm(:, :, s) .* 1e-3) .* sin(phi(:, :, s)) + zg0 + (offcentre * sin((-(s * stpa) * pi / 180) + phase));
            end
        end
    end
phi = squeeze(phi(:, 1, :));
yg = meshgrid(1 : gm2, 1 : gm1);

time = ['The elapsed time for tracking the grits is ', num2str(toc(t)), ' seconds'];
end
