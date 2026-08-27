%% SProfiler
% This function is responsible of calculating the surface roughness

function [raRow, raCol, rmsRow, rmsCol, time] = SRoughness(worked, xinter, xdiv)
    t = tic;
    %xg = min(floor(max(xg, [], 'all') - rds), size(worked, 1));

    % Define the surface height matrix as a 2D array
    xlast = cellfun(@(v) v(end), xinter);
    wheelCenter = mean(xlast, "all");
    [~, wheelCenterIndex] = min(abs(wheelCenter - xdiv));

    surface_heights = worked(:, 1 : wheelCenterIndex);

    % Calculate mean of the surface in row and column directions
    meanRow = mean(surface_heights, 2);
    meanCol = mean(surface_heights, 1);

    % Calculate surface roughness parameters in row direction
    rowDiff = surface_heights - meanRow; % Difference between each element and mean of the column
    raRow = sum(abs(rowDiff), 2) ./ size(surface_heights, 2); % Roughness Average in row direction
    rmsRow = sqrt(sum(rowDiff .^ 2, 2)) ./ size(surface_heights, 2); % Root Mean Square in row direction

    % Calculate surface roughness parameters in column direction
    colDiff = surface_heights - meanCol; % Difference between each element and mean of the row
    raCol = sum(abs(colDiff), 1) ./ size(surface_heights, 1); % Roughness Average in column direction
    rmsCol = sqrt(sum(colDiff .^ 2, 1)) ./ size(surface_heights, 1); % Root Mean Square in column direction

    time = ['The elapsed time for computing the RAs and RMSs is ', num2str(toc(t)), ' seconds'];
end