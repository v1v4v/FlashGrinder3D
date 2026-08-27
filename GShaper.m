function [Heights, gt, gs, time] = GShaper(gm, gs, gsh, gr, xint, yint, gtemp, geq, gmin, gmax, gst)
    t = tic;
    % Grit resolution
    if gst(2) == 1
        disp('Point grit')
        gt = 1;
    elseif gst(2) == 0
        disp('grits are shaped')
        if gsh == "^"
            gt = zeros(gr(1), gr(2));
            vinc = 1 / gr(1);
            hinc = 1 / (floor(gr(2) / 2));
            for r = 1 : 1 : gr(1)
                gt(r, ceil(gr(2) / 2)) = 1 - ((r - 1) * vinc);
            end
            for r = 1 : 1 : gr(1)
                for c = 1 : 1 : floor(gr(2) / 2)
                    gt(r, c) = -(r - 1) * hinc + (c - 0) * vinc;
                    if gt(r, c) < 0
                        gt(r, c) = 0;
                    end
                end
            end
            gt = [gt(:, 1 : 1 : ceil(gr(2) / 2)), flip(gt(:, 1 : 1 : floor(gr(2) / 2)), 2)];
            gt(:, ceil(gr(2) / 2)) = gt(:, ceil(gr(2) / 2)) + hinc;
            gt = flip(gt, 1);
            gt = gt ./ max(gt, [], 'all');
    
    
        elseif gsh == "/"
            gt = zeros(gr(1), gr(2));
            vinc = 1 / gr(1);
            hinc = 1 / gr(2);
            for r = 1 : 1 : gr(1)
                gt(r, gr(2)) = 1 - (r * vinc);
            end
            for r = 1 : 1 : gr(1)-1
                for c = 1 : 1 : gr(2)
                    gt(r, c) = (gr(1) - r + 1) * hinc + (c - (gr(2) / 2)) * vinc;
                    if gt(r, c) < 0
                        gt(r, c) = 0;
                    end
                end
            end
            gt = gt ./ max(gt, [], 'all');
    
        elseif gsh == "\"
            gt = zeros(gr(1), gr(2));
            vinc = 1 / gr(1);
            hinc = 1 / gr(2);
            for r = 1 : 1 : gr(1)
                gt(r, gr(2)) = 1 - ((r - 1) * vinc);
            end
            for r = 1 : 1 : gr(1)
                for c = 1 : 1 : gr(2)
                    gt(r, c) = -(r - 1) * hinc + (c - 0) * vinc;
                    if gt(r, c) < 0
                        gt(r, c) = 0;
                    end
                end
            end
            gt = flip(gt ./ max(gt, [], 'all'), 2);
    
        elseif gsh == "-"
            gt = ones(gr(1), gr(2));
            for r = 1 : 1 : gr(1)
                for c = 1 : 1 : gr(2)
                    if r == 1 || r == gr(1) || c == 1 || c == gr(2)
                        gt(r, c) = .84 * gt(r, c);
                    end
                end
            end
    
            gt = gt ./ max(gt, [], 'all');
    
        elseif gsh == "~"
            syms x y
            geq = str2sym(geq);
            gt = zeros(gr(1), gr(2));
            dx = (xint(2) - xint(1)) / gr(2);
            dy = (yint(2) - yint(1)) / gr(1);
            for xx = 1 : 1 : gr(2)
                xs(xx) = xint(1) + xx * dx;
            end
            for yy = 1 : 1 : gr(1)
                ys(yy) = yint(1) + yy * dy;
            end
            for r = 1 : 1 : gr(1)
                for c = 1 : 1 : gr(2)
                    gt(r, c) = subs(geq, [x, y], [xs(r), ys(c)]);
                end
            end
            gt = gt - min(gt, [], 'all');
            gt = gt ./ max(gt, [], 'all');
    
    
        elseif gsh == "m"
            gt = flip(gtemp, 1);
            %gt = gtemp;
            gt = gt ./ max(gt, [], 'all');
    
        end
        
    end
    % Scales up each elemet in both vertical and horizontal direction based
    % on the dimension of the grit matrix, gr
    Heights = kron(gm, gt);
    scaler = ones(size(gt, 1), size(gt, 2));
    gs = kron(gs, scaler);
    Heights = Heights .* gs;
    time = ['The elapsed time for shaping the grits is ', num2str(toc(t)), ' seconds'];
end