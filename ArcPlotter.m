
function [x, z, theta] = ArcPlotter(pp, ppp, grds, r, dd, s, ad, plott)

    % This function plots an arch passing through the p1 and p2 with a radius of grds
    
    if pp(1) > ppp(1)
        p1 = pp;
        p2 = ppp;
    else
        p1 = ppp;
        p2 = pp;
    end

    % Calculate the distance between the two points
    d = norm(p2 - p1);
    
    % Calculate the center point of the circle
    cx = (p1(1) + p2(1)) / 2 + (p2(2) - p1(2)) / (2 * d) * sqrt((2 * grds) ^ 2 - d ^ 2);
    cz = (p1(2) + p2(2)) / 2 + (p1(1) - p2(1)) / (2 * d) * sqrt((2 * grds) ^ 2 - d ^ 2);
    
    % Calculate the start and end angles of the arc
    theta1 = atan2(p1(2) - cz, p1(1) - cx);
    theta2 = atan2(p2(2) - cz, p2(1) - cx);
    
    % sA = rad2deg(theta1)
    % eA = rad2deg(theta2)
    
    % Create an array of angles for the arc
    % if theta2 < theta1
    %     theta = linspace(theta2, theta1, ad);
    % else
    %     theta = linspace(theta1, theta2, ad);
    % end
    theta = linspace(theta1, theta2, ad);
    % Calculate the x and z coordinates of the arc
    x = cx + grds * cos(theta);
    z = cz + grds * sin(theta);

    %xx = [(pl * ad) + 1 : 1 : (pl + 1) * ad];
    % Plot the arc
    if plott == 'y'
        figure(2)
        hold on
        grid minor
        view(3)
        y = dd * ones(1, ad);
        %work = @(x) 25;
        %fplot(work, [0, 5], 'k', 'LineWidth', 3)
        %plot3(x, y, z, 'LineWidth', 2);
        plot3(x, y, z, LineWidth = 1)
        %axis equal;
        text(mean(x), mean(y), mean(z), strcat(int2str(r), 's', int2str(s)))
        grid on
        drawnow
    end
    theta = rad2deg(theta);
end