clc
clearvars -except PTh_range SA_range vars_file
t = tic;
disp('Scripta manent')
%% Stage 1, Geometrical and process parameters
% Open the file for reading
InputFileName = 'vars2.txt';
fileID = fopen(InputFileName, 'r');

% Read the contents of the file
contents = textscan(fileID, '%s%s', 'Delimiter', '=', 'Whitespace', '\n');

% Extract the variable names and values
variableNames = contents{1};
variableValues = contents{2};
% Evaluates the variables
for i = 1 : numel(variableNames)
    evalc([variableNames{i} ' = ' num2str(variableValues{i})]);
end
fclose(fileID);
% save([MyDir, '\testInputs', num2str(TestNo), '.txt'], "-ascii")
clear contents variableNames variableValues valueStr fileID ans
symm = symmetry{symm};
%% Stage 2, Material parameters

% Open the text file for reading WORKPIECE properties
mtrl = mtrlz{mtrl};
fileID = fopen('mtrlz.txt', 'r');
% Read the file line by line
while ~feof(fileID)
    line = fgetl(fileID);
    target = append('~', mtrl);
    % Check if the line starts with the delimiter
    if strncmp(line, target, 9)
        % Initialize a structure to store the material properties
        material = struct('Material', mtrl, 'A', 0, 'B', 0, 'n', 0, 'm', 0, 'C', 0, 'Tm', 0, 'rho', 0, 'k', 0, 'cp', 0);
        % Read the data until the next delimiter or end of file
        while ~feof(fileID)
            line = fgetl(fileID);
            % Removes the whitespaces
            line = line(~isspace(line));                
            % Check if the line starts with a property identifier followed by '='
            if contains(line, '=')
                % Split the line into property name and value
                contents = textscan(line, '%s%s', 'Delimiter', '=', 'Whitespace', '\n');    
                material.(string(contents{1})) = eval(cell2mat(contents{2}));                    
            elseif strncmp(line, '~', 1)
                break;
            end
        end
    end
end
% Close the file
fclose(fileID);
clear contents line target

% Open the text file for reading COOLANT properties
clnt = clntz{clnt};
fileID = fopen('clntz.txt', 'r');
% Read the file line by line
while ~feof(fileID)
    line = fgetl(fileID);
    target = append('~', clnt);
    % Check if the line starts with the delimiter
    if strncmp(line, target, 7)
        % Initialize a structure to store the material properties
        coolant = struct('Coolant', clnt, 'k', 0, 'cp', 0, 'rho', 0, 'mu', 0, 'nu', 0, 'tboil', 0);
        % Read the data until the next delimiter or end of file
        while ~feof(fileID)
            line = fgetl(fileID);
            % Removes the whitespaces
            line = line(~isspace(line));                
            % Check if the line starts with a property identifier followed by '='
            if contains(line, '=')
                % Split the line into property name and value
                contents = textscan(line, '%s%s', 'Delimiter', '=', 'Whitespace', '\n');    
                coolant.(string(contents{1})) = eval(cell2mat(contents{2}));                    
            elseif strncmp(line, '~', 1)
                break;
            end
        end
    end
end
% Close the file
fclose(fileID);
clear contents line target

% Open the text file for reading GRIT properties
wheel = wheelz{wheel};
grt = grtz{grt};
fileID = fopen('grtz.txt', 'r');
% Read the file line by line
while ~feof(fileID)
    line = fgetl(fileID);
    target = append('~', grt);
    % Check if the line starts with the delimiter
    if strncmp(line, target, 7)
        % Initialize a structure to store the material properties
        grit = struct('Grit', grt, 'k', 0, 'cp', 0, 'rho', 0, 'mu', 0, 'nu', 0, 'tboil', 0);
        % Read the data until the next delimiter or end of file
        while ~feof(fileID)
            line = fgetl(fileID);
            % Removes the whitespaces
            line = line(~isspace(line));                
            % Check if the line starts with a property identifier followed by '='
            if contains(line, '=')
                % Split the line into property name and value
                contents = textscan(line, '%s%s', 'Delimiter', '=', 'Whitespace', '\n');    
                grit.(string(contents{1})) = eval(cell2mat(contents{2}));                    
            elseif strncmp(line, '~', 1)
                break;
            end
        end
    end
end
% Close the file
fclose(fileID);
clear contents line target

%% Saving the inputs
outputFolder = fullfile(cd, MyDir);
if ~exist(outputFolder, 'dir')
    mkdir(MyDir);
end
destFile = fullfile(outputFolder, ['testInputs', num2str(TestNo), '.txt']);
copyfile(fullfile(cd, InputFileName), destFile)
% save([MyDir, '\testInputs', num2str(TestNo), '.mat'])
close all
%% Stage 3, Fully-Coupled simulator
% OneGrit creates the patterened wheel and computes the trajectory of the grits, the ground work surface, forces, and temperature
OneGrit

%% Stage 4, Storing the results
% computes the overall elapsed time. Works only if OneGrit is already called.
deltat{end + 1} = ['The overall elapsed time is ', num2str(toc(t)), ' seconds'];

% pcode("OneGrit.m", "GPatterner.m", "GSpecifier.m", "GShaper.m", "GParameters.m", "GTracker4.m", "WIntersections.m", "WSurface2.m", "ArcPlotter.m", "GForce1.m", "GForce2.m", "GTempFD.m", "GTempFD3D3.m", "SRoughness.m")
