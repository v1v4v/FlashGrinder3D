

function [xinter, zinter, time] = WIntersections(xg, zg, gm, pg, wh, gmax, fm)

    t = tic;
    xinter = cell(size(xg, 1), size(xg, 2));
    zinter = cell(size(xg, 1), size(xg, 2));
    xinterTot = cell(size(xg, 1), size(xg, 2));
    zinterTot = cell(size(xg, 1), size(xg, 2));

    [xinter{:}] = deal(zeros(size(xg, 3), 1));
    [zinter{:}] = deal(zeros(size(zg, 3), 1));

    limit = wh + (gmax * 1e-3);
    for c = 1 : size(xg, 2)
        for rr = 1 : 1 : length(pg{c})
            r = pg{c}(rr);

            xinterTot{r, c} = squeeze(xg(r, c, :));
            zinterTot{r, c} = squeeze(zg(r, c, :));

            xinter{r, c} = xinterTot{r, c};
            zinter{r, c} = zinterTot{r, c};

            index = find(zinterTot{r, c} < limit);

            xinter{r, c}(index) = xinterTot{r, c}(index);
            zinter{r, c}(index) = zinterTot{r, c}(index);
                        
        end
    end

    % for c = 1 : 1 : size(xg, 2)
    %     for rr = 1 : 1 : length(pg{c}) - 1
    %         r1 = pg{c}(rr);
    %         r2 = pg{c}(rr + 1);
    %         [xints{r, c}, zints{r, c}] = intersections(xinter{r1, c}, zinter{r1, c}, xinter{r2, c}, zinter{r2, c});
    %     end
    % end
    time = ['The elapsed time for finding the intersections of trajectories is ', num2str(toc(t)), ' seconds'];
end