%% GPatterner
% This module creates the full wheel pattern
% first a template is created and then based on the desired pattern, it is
% replicated horizontally and vertically


function [Elements, fullpattern, remRatio, time] = GPatterner(WheelDiameter, WheelWidth, StripesNumber, StripesAngle, Fn1, Fn2, PT, eq, int, PatternDir, pth, symm, ps, gd, wheel, gst)
    %  Discretization of the wheel peripheral surface
    %  Between each particle is divided to 'Fineness' elements
    %     - WheelDiameter: Wheel diameter [cm]
    %     - WheelWidth: Wheel width [cm]
    %     - StripesNumber: Number of pattern repetitions in the width of wheel [#]
    %     - StripesAngle: Angle of the pattern for linear patterns, i.e., '/', 'v', '^', 'x' [degree]
    %     - Fn1: Number of elements in the width of wheel for each patten repetition  [#]
    %     - Fn2: Total number of elements in the periphery of the wheel [#]
    %     - PT: Pattern type
    %     - eq: Pattern equation for '~' type
    %     - int: Interval for the equation of '~' pattern type 
    %     - PatternDir: Pattern replication direction, 'hor' for replication in the width, and 'ver' for replication in the peripheral direction
    %     - pth: Pattern thickness [number of elements]
    %     - symm: Symmetrical pattern creation for '~' pattern type
    %     - ps: Spacing between each full pattern row [number of elements]
    %     - gd: Grit density on the wheel surface [%]
    %     - wheel: Wheel to be "generated", "loaded" or "loaded and patterned"
    %     - gst: General simulation tweaks
    %
    %     - Elements: Final state of dicretized wheel either with pattern or not
    %     - fullpattern: The final version of the wheel without the considering the grit density
    %     - remRatio: Remaining ratio of the grits in [mm^2]
    %     - time: Elapsed time


    t = tic;
    if pth > 0
        disp("Wheel having " + PT + " pattern with " + num2str(pth) + " elements thickness.")
    else
        disp("Full wheel.")
    end
    ElementNumberWidth = (StripesNumber) * Fn1;
    
    
    qqq = ceil(Fn1 * tand(StripesAngle));
    

    Elements = zeros(Fn2, ElementNumberWidth);
    ElementsTemplate = zeros(qqq, Fn1);
    ElementsTemplateAspectRatio = ceil(Fn1 / qqq);
    Degrees = zeros(qqq, Fn1);
    DegreesDifference = ones(qqq, Fn1) * 90;
    Leastd = cell(qqq, Fn1);
    if pth > 0
        if PT == '/' || PT == '\' || PT == 'v' || PT == '^' || PT == 'x'
            disp(PT)                
            for q = 1 : qqq
                if q == 1
                    Leastd{q, Fn1} = StripesAngle;
                    ElementsTemplate(q, Fn1 - ElementsTemplateAspectRatio + 1 : Fn1) = 1;
                else
                    for w = 1 : Fn1
                        Degrees(q, w) = atand((q - 1) / (Fn1 - w ));
                        if abs(StripesAngle - Degrees(q, w)) < DegreesDifference(q, w)
                            DegreesDifference(q, w) = abs(StripesAngle - Degrees(q, w));
                        end
                    end
                    [mn, ind] = sort(DegreesDifference(q, :));
                    Leastd{q, ind(1)} = mn(1 : ElementsTemplateAspectRatio);
                    ElementsTemplate(q, ind(1 : ElementsTemplateAspectRatio)) = 1;
                end
            end
        elseif PT == "~"
            disp(eq)
            eq = str2sym(eq);
            xs = linspace(int(1), int(2), Fn1);
            ys = eval(subs(eq, xs));
            ys = ys ./ max(ys);
            ys = ys .* 10;
            yspan = ceil(max(ys) - min(ys));
            ElementsTemplate = zeros(yspan, Fn1);        
            for c = 1 : 1 : Fn1
                if min(ys) < 0
                    ElementsTemplate(round(ys(c) + abs(min(ys)) + 1), c) = 1;
                else
                    ElementsTemplate(floor(ys(c) - abs(min(ys)) + 1), c) = 1;
                end
            end
    
            switch lower(symm)
                case 'verticalnhorizontal'
                    ElementsTemplate = [ElementsTemplate; flip(ElementsTemplate, 1)];
                    ElementsTemplate = [ElementsTemplate, flip(ElementsTemplate, 2)];
                case 'vertical'
                    ElementsTemplate = [ElementsTemplate; flip(ElementsTemplate, 1)];
                case 'horizontal'
                    ElementsTemplate = [ElementsTemplate flip(ElementsTemplate, 2)];
                otherwise
                    disp('Asymmetric pattern.');
            end
    
        end

        if size(ElementsTemplate, 2) > Fn1
            ElementsTemplate = double(imbinarize(imresize(ElementsTemplate, Fn1 / size(ElementsTemplate, 2))));
        end
        if any(strcmp(PT, {'x'}))
            ElementsTemplate = imbinarize(ElementsTemplate + flip(ElementsTemplate, 2));
        end
        if any(strcmp(PT, {'v', '^'}))
            ElementsTemplate = double(imbinarize(imresize(ElementsTemplate, .5)));
        end

        
        PeripheralMultiplication = ceil(Fn2 / size(ElementsTemplate, 1));
        elte = ElementsTemplate;
        dis = pth + ps;
        if dis < size(ElementsTemplate, 1)
            elte2 = [elte; zeros(Fn2, size(elte, 2))];
            elte3 = elte2;
            unitt = 1;
            displus = dis + 1;
            while unitt * displus < 2 * size(elte, 1)
                elte3 = imbinarize(elte3 + circshift(elte2, unitt * displus, 1));
                unitt = unitt + 1;
            end
            ElementsTemplate = elte3(size(elte, 1) + 1 : 2 * size(elte, 1), :);
        elseif dis > size(ElementsTemplate, 1)
            tobeadd = abs(dis - size(ElementsTemplate, 1));
            ElementsTemplate = [elte; zeros(tobeadd + 1, size(elte, 2))];
            % ElementsTemplate = imbinarize(elte2 + circshift(elte2, dis + 1, 1)+ circshift(elte2, 2 * (dis + 1), 1));
            % ElementsTemplate = ElementsTemplate(1 : 2 * (dis + 1) + size(elte, 1) + abs(dis - size(elte, 1)), :);
        end

        if gst(1) == 1
            spacing = ones(0, size(ElementsTemplate, 2));
        elseif gst(1) == 0
            spacing = zeros(0, size(ElementsTemplate, 2));
        end
        
    end

    % if ~strcmp(PT, '~') && ~strcmp(PT, 'x')
    %     s2e = strel('line', ceil(pth * sind(StripesAngle)), (StripesAngle + 90));
    %     s3e = strel('square', 2);
    %     s2e = strel(imdilate(s2e.Neighborhood, s3e));
    % elseif strcmp(PT, '~') || strcmp(PT, 'x')
    %     s2e = strel('disk', floor(pth / 2));
    % end

    if ~strcmp(PT, '~') && ~strcmp(PT, 'x')
        s2e = strel('line', ceil(pth / sind(90 - StripesAngle)), 90);
        % s3e = strel('square', 2);
        % s2e = strel(imdilate(s2e.Neighborhood, s3e));
        % s2e = strel('disk', floor(pth / 2));
    elseif strcmp(PT, '~') || strcmp(PT, 'x')
        s2e = strel('line', ceil(pth / 2), 90);
    end

    if ~strcmp(PT, '~') && ~strcmp(PT, 'x')
        ElementsTemplate = [zeros(ceil((pth / 2 / sind(90 - StripesAngle))), size(ElementsTemplate, 2)); ElementsTemplate; zeros(ceil((pth / 2 / sind(90 - StripesAngle))), size(ElementsTemplate, 2))];
    elseif strcmp(PT, '~') || strcmp(PT, 'x')
        ElementsTemplate = [zeros(round(pth / 2), size(ElementsTemplate, 2)); ElementsTemplate; zeros(round(pth / 2), size(ElementsTemplate, 2))];
    end

    ElementsTemplate = imdilate(ElementsTemplate, s2e);

    Elements = ElementsTemplate;
    % Elements = imdilate(Elements, s2e);
    if size(Elements, 2) > (Fn1 * StripesNumber) || size(Elements, 1) > Fn2
        if size(Elements, 2) > (Fn1 * StripesNumber)
            trim2 = size(Elements, 2) - (Fn1 * StripesNumber);
        else
            trim2 = 0;
        end
        if size(Elements, 1) > Fn2
            trim1 = size(Elements, 1) - Fn2;
        else
            trim1 = 0;
        end
        Elements = Elements(ceil(trim1 / 2) + 1 : end - ceil(trim1 / 2), ceil(trim2 / 2) + 1 : end - floor(trim2 / 2));
    end
    Elements = double(imbinarize(imresize(Elements, [Fn2 (Fn1 * StripesNumber)])));
    

    if gst(1) == 1
        if pth < 1
            Elements(Elements == 1) = 0;
        end
        Elements(Elements == 1) = 2;
        Elements(Elements == 0) = 1;
        Elements(Elements == 2) = 0;
    end
    fullpattern = Elements;
    if strcmp(wheel, "create")
        for r = 1 : size(Elements, 1)
            idx = find(Elements(r, :) ~= 0);
            rem = ceil(gd * length(idx));
            fallengrits = randsample(idx, length(idx) - rem);
            Elements(r, fallengrits) = 0;
        end
    end
    % Ratio of remaining area
    remRatio = numel(nonzeros(fullpattern)) / numel(fullpattern);
    if pth > 0
        disp(['Ratio of remaining area after laser ablation is = ', num2str(remRatio)])
    else
        disp(['Ratio of remaining area after laser ablation is = ', num2str(1)])
    end
    %writematrix(Elements, "Pattern.xlsx",'Sheet',1);
    time = ['The elapsed time for patterning is ', num2str(toc(t)), ' seconds'];
end