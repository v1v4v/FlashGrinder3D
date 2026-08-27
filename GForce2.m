%% GForce2

function [fx, fy, fz, fxi1, fxi12, fxi13, fzi1, fzi12, fzi13, percentplough, percentcut, percentdmz, ft, fn, thm, time] = GForce2(ad, hcr, chth, inc, p0, px, tau, signz, l1, l2, l3, mu, kt, kn, kr, gs, rdoc, theta, steptime)
    t = tic;
    % i means infinitesimal, or parfor each single grit
    % x as the last charachter in names, means extended, from cell array to 3D matrix, i.e. fycx, chthx, etc.
    % 1, 2, and 3 refer to the three flank face lengths


    ffx1 = zeros(size(hcr));
    ffz1 = zeros(size(hcr));
    ffx2 = zeros(size(hcr));
    ffz2 = zeros(size(hcr));
    ffx3 = zeros(size(hcr));
    ffz3 = zeros(size(hcr));

    fnx1 = zeros(size(hcr));
    fnz1 = zeros(size(hcr));
    fnx2 = zeros(size(hcr));
    fnz2 = zeros(size(hcr));
    fnx3 = zeros(size(hcr));
    fnz3 = zeros(size(hcr));

    
    ftc = cell(size(hcr));
    fnc = cell(size(hcr));
    frc = cell(size(hcr));
    
    fxc = cell(size(hcr));
    fyc = cell(size(hcr));
    fzc = cell(size(hcr));

    ftdmz = cell(size(hcr));
    fadmz = cell(size(hcr));
    fndmz = cell(size(hcr));

    fxdmz = cell(size(hcr));
    fydmz = cell(size(hcr));
    fzdmz = cell(size(hcr));


    dt = cell(size(hcr));

    chthx = zeros(size(hcr, 1), size(hcr, 2), ad);
    thetax = zeros(size(hcr, 1), size(hcr, 2), ad);

    fxcx = zeros(size(hcr, 1), size(hcr, 2), ad);
    fycx = zeros(size(hcr, 1), size(hcr, 2), ad);
    fzcx = zeros(size(hcr, 1), size(hcr, 2), ad);

    fxdmzx = zeros(size(hcr, 1), size(hcr, 2), ad);
    fydmzx = zeros(size(hcr, 1), size(hcr, 2), ad);
    fzdmzx = zeros(size(hcr, 1), size(hcr, 2), ad);


    % Ploughing force construction
    for d = 1 : size(gs, 2)
        for r = inc{1, d}
            % Friction
            zone31 = linspace(0, l1(r, d), 100);
            ffx1(r, d) = trapz(zone31, heaviside(-signz{r, d}) .* (mu(r, d) .* px{r, d} .* sind(90 + (zone31 ./ (gs(r, d, 5) * 1e-3)) - gs(r, d, 8)))) + trapz(zone31, heaviside(signz{r, d}) .* (tau(r, d) .* gs(r, d, 4) .* sind(90 + (zone31 ./ (gs(r, d, 5) * 1e-3)) - gs(r, d, 8))));
            ffz1(r, d) = trapz(zone31, heaviside(-signz{r, d}) .* (mu(r, d) .* px{r, d} .* cosd(90 + (zone31 ./ (gs(r, d, 5) * 1e-3)) - gs(r, d, 8)))) + trapz(zone31, heaviside(signz{r, d}) .* (tau(r, d) .* gs(r, d, 4) .* cosd(90 + (zone31 ./ (gs(r, d, 5) * 1e-3)) - gs(r, d, 8))));
            
            zone32 = linspace(l1(r, d), (l1(r, d) + l2(r, d)), 100);
            ffx2(r, d) = trapz(zone32, heaviside(-signz{r, d}) .* (mu(r, d) .* px{r, d} .* cosd((zone32 ./ (gs(r, d, 5) * 1e-3)) - gs(r, d, 8)))) + trapz(zone32, heaviside(signz{r, d}) .* (tau(r, d) .* gs(r, d, 4) .* cosd((zone32 ./ (gs(r, d, 5) * 1e-3)) - gs(r, d, 8)))) ;
            ffz2(r, d) = trapz(zone32, heaviside(-signz{r, d}) .* (mu(r, d) .* px{r, d} .* sind((zone32 ./ (gs(r, d, 5) * 1e-3)) - gs(r, d, 8)))) - trapz(zone32, heaviside(signz{r, d}) .* (tau(r, d) .* gs(r, d, 4) .* sind((zone32 ./ (gs(r, d, 5) * 1e-3)) - gs(r, d, 8)))) ;
            
            zone33 = linspace((l1(r, d) + l2(r, d)), (l1(r, d) + l2(r, d) + l3(r, d)), 100);
            ffx3(r, d) = trapz(zone33, heaviside(-signz{r, d}) .* (mu(r, d) .* px{r, d} .* cosd(gs(r, d, 3)))) + trapz(zone33, heaviside(signz{r, d}) .* tau(r, d) * gs(r, d, 4) .* cosd(gs(r, d, 3)));
            ffz3(r, d) = - trapz(zone33, heaviside(-signz{r, d}) .* (mu(r, d) .* px{r, d} .* sind(gs(r, d, 3)))) - trapz(zone33, heaviside(signz{r, d}) .* tau(r, d) * gs(r, d, 4) .* sind(gs(r, d, 3)));
            
            % Normal
            fnx1(r, d) = trapz(zone31, (gs(r, d, 4) .* px{r, d} .* cosd(90 + (zone31 ./ (gs(r, d, 5) * 1e-3)) - gs(r, d, 8))));
            fnz1(r, d) = - trapz(zone31, (gs(r, d, 4) .* px{r, d} .* sind(90 + (zone31 ./ (gs(r, d, 5) * 1e-3)) - gs(r, d, 8))));
        
            fnx2(r, d) = trapz(zone32, (gs(r, d, 4) .* px{r, d} .* sind((zone32 ./ (gs(r, d, 5) * 1e-3)) - gs(r, d, 8))));
            fnz2(r, d) = - trapz(zone32, (gs(r, d, 4) .* px{r, d} .* cosd((zone32 ./ (gs(r, d, 5) * 1e-3)) - gs(r, d, 8))));
        
            fnx3(r, d) = - trapz(zone33, (gs(r, d, 4) .* px{r, d} .* sind(gs(r, d, 3))));
            fnz3(r, d) = - trapz(zone33, (gs(r, d, 4) .* px{r, d} .* cosd(gs(r, d, 3))));
        end
    end

    ffx = ffx1 + ffx2 + ffx3;
    ffz = ffz1 + ffz2 + ffz3;

    fnx = fnx1 + fnx2 + fnx3;
    fnz = fnz1 + fnz2 + fnz3;

    % Forces on the tool
    fxp = (ffx + fnx);
    fyp = 0;
    fzp = (ffz + fnz);

    % Cutting force construction
    
    for d = 1 : 1 : size(hcr, 2)
        for r = inc{1, d}
            chthx(r, d, 1 : ad) = chth{r, d};
            thetax(r, d, 1 : ad) = theta{r, d};
        end
    end
    
    normalizedchth = cellfun(@(v) v / max(v), chth, 'UniformOutput', false);

    for d = 1 : 1 : size(gs, 2)
        for r = inc{:, d}
            ftc{r, d} = (1 / 2) .* gs(r, d, 4) .* kt(r, d) .* chth{r, d} .* normalizedchth{r, d};
            fnc{r, d} = (1 / 2) .* gs(r, d, 4) .* kn(r, d) .* chth{r, d} .* normalizedchth{r, d};
            frc{r, d} = (1 / 2) .* gs(r, d, 4) .* kr(r, d) .* chth{r, d} .* normalizedchth{r, d};
            %ftcave(r, d) = mean(ftc{r, d});

            fxc{r, d} = ftc{r, d} .* cosd(abs(90 - theta{r, d})) - fnc{r, d} .* cosd(abs(theta{r, d}));
            fyc{r, d} = 0;
            fzc{r, d} = - ftc{r, d} .* sind(abs(90 - theta{r, d})) - fnc{r, d} .* sind(abs(theta{r, d}));
        end
    end


    for d = 1 : 1 : size(hcr, 2)
        for r = inc{1, d}
            fxcx(r, d, 1 : ad) = fxc{r, d};
            fycx(r, d, 1 : ad) = fyc{r, d};
            fzcx(r, d, 1 : ad) = fzc{r, d};
        end
    end
    
    % DMZ force construction
    % Uncut chip thickness is set to be RDOC, but in reality it is different for each grit and must be modified! Has been modified!!

    for d = 1 : 1 : size(hcr, 2)
        for r = inc{1, d}
            dt{r, d} = chth{r, d} .* deg2rad(abs(gs(r, d, 1)));
            % dt{r, d} = chth{r, d} .* (tand(abs(gs(r, d, 1))) .^ 1.64);
            ftdmz{r, d} = tau(r, d) .* dt{r, d} .* (gs(r, d, 4)) .* normalizedchth{r, d};
            fadmz{r, d} = 0;
            fndmz{r, d} = p0(r, d) .* dt{r, d} .* (gs(r, d, 4)) .* normalizedchth{r, d};

            fxdmzx(r, d, 1 : ad) = (ftdmz{r, d} .* cosd(abs(90 - theta{r, d}))) - (fndmz{r, d} .* cosd(abs(theta{r, d})));
            fydmzx(r, d, 1 : ad) = fadmz{r, d};
            fzdmzx(r, d, 1 : ad) = - (ftdmz{r, d} .* sind(abs(90 - theta{r, d}))) - (fndmz{r, d} .* sind(abs(theta{r, d})));
        end
    end


    % for d = 1 : 1 : size(hcr, 2)
    %     for r = inc{1, d}
    %         fxdmzx(r, d, 1 : ad) = ftdmz{r, d};
    %         fydmzx(r, d, 1 : ad) = fadmz{r, d};
    %         fzdmzx(r, d, 1 : ad) = fndmz{r, d};
    %     end
    % end
    
    index = ~cellfun(@isempty, theta);
    thetamean = cellfun(@mean, theta);
    thm = 90 - mean(thetamean(index));
    thetax = 90 - thetax;

    comparison1 = chthx > hcr;
    comparison2 = chthx < hcr & chthx > 0;
    sum1 = sum(comparison1, 3);
    sum2 = sum(comparison2, 3);

    fxi1 = (trapz(comparison2 .* fxp, 3) ./ sum2);
    fxi12 = (trapz(comparison1 .* (fxcx), 3) ./ sum1);
    fxi13 = (trapz(comparison1 .* (fxdmzx), 3) ./ sum1);
    fxi2 = (trapz(comparison1 .* (fxcx + fxdmzx), 3) ./ sum1);
    fxit = (trapz(comparison1 .* ((fxcx + fxdmzx + fxp) .* cosd(thetax)), 3) ./ sum1) + (trapz(comparison2 .* fxp .* cosd(thetax), 3) ./ sum2);
    fxin = (trapz(comparison1 .* ((fxcx + fxdmzx + fxp) .* sind(thetax)), 3) ./ sum1) + (trapz(comparison2 .* fxp .* sind(thetax), 3) ./ sum2);

    fxi1(isnan(fxi1)) = 0;
    fxi12(isnan(fxi12)) = 0;
    fxi13(isnan(fxi13)) = 0;
    fxi2(isnan(fxi2)) = 0;
    fxit(isnan(fxit)) = 0;
    fxin(isnan(fxin)) = 0;
    
    fyi1 = (trapz(comparison2 .* fyp, 3) ./ sum2);
    fyi2 = (trapz(comparison1 .* (fycx + fydmzx), 3) ./ sum1);
    fyi1(isnan(fyi1)) = 0;
    fyi2(isnan(fyi2)) = 0;

    fzi1 = (trapz(comparison2 .* fzp, 3) ./ sum2);
    fzi12 = (trapz(comparison1 .* (fzcx), 3) ./ sum1);
    fzi13 = (trapz(comparison1 .* (fzdmzx), 3) ./ sum1);
    fzi2 = (trapz(comparison1 .* (fzcx + fzdmzx), 3) ./ sum1);
    fzit = (trapz(comparison1 .* ((fzcx + fzdmzx + fzp) .* sind(thetax)), 3) ./ sum1) + (trapz(comparison2 .* fzp .* sind(thetax), 3) ./ sum2);
    fzin = (trapz(comparison1 .* ((fzcx + fzdmzx + fzp) .* cosd(thetax)), 3) ./ sum1) + (trapz(comparison2 .* fzp .* cosd(thetax), 3) ./ sum2);

    fzi1(isnan(fzi1)) = 0;
    fzi12(isnan(fzi12)) = 0;
    fzi13(isnan(fzi13)) = 0;
    fzi2(isnan(fzi2)) = 0;
    fzit(isnan(fzit)) = 0;
    fzin(isnan(fzin)) = 0;

    fxi = fxi1 + fxi12 + fxi13;
    fyi = fyi1 + fyi2;
    fzi = fzi1 + fzi12 + fzi13;

    ft = fxit + fzit;
    fn = fxin + fzin;
    
    % Total Ploughing Force
    fxi1 =  sum(fxi1, 'all');
    fzi1 =  sum(fzi1, 'all');

    % Total cutting Force
    fxi12 =  sum(fxi12, 'all');
    fzi12 =  sum(fzi12, 'all');

    % Total DMZ Force
    fxi13 =  sum(fxi13, 'all');
    fzi13 =  sum(fzi13, 'all');

    % Total force
    fx = sum(fxi, 'all');
    fy = sum(fyi, 'all');
    fz = sum(fzi, 'all');

    % Percentage of components
    percentplough(1, 1) = abs(fxi1 / fx);
    percentplough(1, 3) = abs(fzi1 / fz);
    
    percentcut(1, 1) = abs(fxi12 / fx);
    percentcut(1, 3) = abs(fzi12 / fz);

    percentdmz(1, 1) = abs(fxi13 / fx);
    percentdmz(1, 3) = abs(fzi13 / fz);

    ft = sum(ft, 'all');
    fn = sum(fn, 'all');
    if abs(ft) < .1
        ft = sign(ft) * .1;
    end

    time = ['The elapsed time for computing the force components is ', num2str(toc(t)), ' seconds'];
end