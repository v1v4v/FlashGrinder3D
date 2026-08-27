function [params] = GTempFD(timstart, timend, subtimeres, wcl, workdimz, heightMap, pretemp, heatPattern, AreaL, div, xg, f, tamb, tfluid, material, fdmethod, coolant, grit, rdoc, wd, ww, wr, gm, vc, epsilon, bcs, deltmax, tch, fdm)


    % Finite Difference Solver parfor any transient 2 or 3D property (Forward-Time, Central-Space, MacCormack, Crank-Nicolson)
    
    % End of the time-step
    timend;
    % Start of the time-step
    timstart;

    % Force per unit width of the wheel [N/mm]
    f = abs(f);
    tamb = tamb;
    tfluid = tfluid;
    % tch = max(tempmax, 150);
    if abs(f) < .01
        f = .01;
    else
        f = abs(f);
    end
    x = mean(xg(:, :), 'all');
    Lc = sqrt(rdoc * wd);
    %grit
    alphag = grit.k / (grit.rho * grit.cp);
    alphac = coolant.k / (coolant.rho * coolant.cp);
    alpham = material.k / (material.rho * material.cp);

    phi = sqrt((alphag .* Lc) ./ (((wr + (gm .* 1e-3)) .^ 2) .* fdm));
    F = 1 - exp(- phi ./ 1.2);
    Rws = (1 + ((grit.k) ./ (sqrt(material.k .* material.rho .* material.cp .* (wr + (gm * 1e-3)) .* vc) .* F))) .^ -1;

    % coolant.mu is viscousity of the coolant fluid, called eta in the paper
    h = (4 / 9) .* ((coolant.k) .^ (2 / 3)) .* ((coolant.rho) .^ (1 / 2)) .* ((coolant.cp) .^ (1 / 3)) .* ((coolant.mu) .^ (- 1 / 6)) .* sqrt(vc ./ Lc);
    hmean = mean(h, "all") .* 1e6;

    
    InfluxArea = workdimz(1) * Lc;
    % Heat out flux from the head of the workpiece
    OutfluxArea1 = workdimz(1) * workdimz(2);
    % Heat out flux from the sides of the workpiece
    OutfluxArea2 = workdimz(2) * workdimz(3);
    % Heat out flux from the top of the workpiece till wheel
    OutfluxArea3 = workdimz(1) * wcl;
    
    qys = linspace(0, workdimz(3), div);
    
    [~, cnt(1)] = min(abs(wcl - qys));
    [~, cnt(2)] = min(abs((wcl + Lc) - qys));

    % InfluxElements = abs(cnt(2) - cnt(1)) * workdimz(1);
    % qtot = qch + qw + qwh + qf;
    qch = (rdoc * fdm * material.rho * material.cp * tch) / Lc;
    qtot = (f * mean(vc, 'all') * 1e-3);
    qf1 = -mean(h, 'all') * abs(tch - tfluid);
    qw = Rws * (qtot - qch) - qf1;

    eps = max(qw, [], 'all') / qtot;
    
    divi = size(pretemp, 3);
    delt = (timend - timstart) / subtimeres;
    delx = workdimz(3) / div;
    dely = workdimz(1) / size(gm, 2);
    delz = workdimz(2) / divi;

    penet = 1;
    % convo1 = sin(linspace(0, pi, cnt(2) - cnt(1) + 1)); 
    convo1 = 1;
    % convo2 = repmat(cos(linspace(0, pi / 2, round(penet * cnt(1)))), size(workfield, 1), 1);
    convo2 = 1;
    %  material.rho * material.cp *
    % * AreaL
    qin = (epsilon * qtot) / ((InfluxArea * 1e-6));
    qwmax = max(qw, [], 'all') / ((InfluxArea * 1e-6));
    % qin(2 : end - 1, cnt(1) : cnt(2)) = (epsilon * qtot) / (material.rho * material.cp * InfluxArea * AreaL);
    % qin(2 : end - 1, cnt(1) : cnt(2)) = qin(2 : end - 1, cnt(1) : cnt(2)) .* convo1;

    
    % qin(2 : end - 1, 1 : round(penet * cnt(1))) = qin(2 : end - 1, 1 : round(penet * cnt(1))) + qf2(2 : end - 1, :) .* convo2;
    
    % bcval = max(qwmax, [], "all");
    % bcs.z1 = struct('type', 'Neumann', 'value', -bcval);

    % Parameters
    nx = div;
    ny = size(gm, 2);
    nz = divi;
    numSteps = subtimeres;
    kappa = alpham;

    %% 1. Define Simulation Parameters
    params = struct();
    
    % --- Method
    params.method = fdmethod;
    % --- Physical Dimensions (in MILLIMETERS) ---
    params.DIMS.Ly = workdimz(3) * 1e-3; 
    params.DIMS.Lx = workdimz(1) * 1e-3; 
    params.DIMS.Lz = workdimz(2) * 1e-3;
    
    % --- Grid Resolution ---
    params.GRID.Ny = ny;
    params.GRID.Nx = nx;
    params.GRID.Nz = nz;
    
    % --- Material Properties (in standard SI units) ---
    params.MATERIAL.alpha = alpham * 1e-6; % Thermal diffusivity (m^2/s)
    params.MATERIAL.k     = material.k * 1e3;    % Thermal conductivity (W/mK)
    
    % --- Time Control ---
    params.TIME.t_final = timend - timstart; % Total simulation time in seconds
    params.TIME.dt      = delt;  % Let the solver calculate a stable time step
    
    % --- Initial Condition (in CELSIUS) ---
    params.T_initial = pretemp;
    
    %% 2. Define Boundary Conditions (Temperatures in CELSIUS)
    
    % Assumptions:
    %   - bcs is a struct with fields: x0,x1,y0,y1,z0,z1
    %   - Each bcs.(name) has at least: .type and .value
    %   - Types can be: "Dirichlet", "Neumann", "Convective", "Mixed" (any case)
    
    bcNames = {'x0','x1','y0','y1','z0','z1'};
    
    if ~exist('params','var') || ~isstruct(params)
        params = struct();
    end
    if ~isfield(params,'BC') || ~isstruct(params.BC)
        params.BC = struct();
    end
    
    for k = 1:numel(bcNames)
        name = bcNames{k};
    
        if ~isfield(bcs, name)
            error('bcs is missing field "%s".', name);
        end
    
        bc = bcs.(name);
    
        if ~isfield(bc,'type')
            error('bcs.%s is missing the field "type".', name);
        end
    
        % Normalize type to lower-case char
        bcType = lower(string(bc.type));
        params.BC.(name).type = char(bcType);
    
        switch bcType
            case "dirichlet"
                % Example: params.BC.z0.type='dirichlet'; params.BC.z0.value=37.0;
                if ~isfield(bc,'value')
                    error('bcs.%s (Dirichlet) must contain .value (scalar).', name);
                end
                params.BC.(name).value = bc.value;
    
            case "neumann"
                % Example: params.BC.y0.type='neumann'; params.BC.y0.value=0.0;
                if ~isfield(bc,'value')
                    error('bcs.%s (Neumann) must contain .value (scalar).', name);
                end
                params.BC.(name).value = bc.value;
    
            case "convective"
                % Example:
                % params.BC.x0.type='convective';
                % params.BC.x0.value.h = hmean;
                % params.BC.x0.value.T_inf = 22.0;
                %
                % Accept either bc.value.h / bc.value.T_inf OR bc.h / bc.T_inf
                if isfield(bc,'value') && isstruct(bc.value)
                    v = bc.value;
                else
                    v = struct();
                    if isfield(bc,'h'),     v.h     = bc.h;     end
                    if isfield(bc,'T_inf'), v.T_inf = bc.T_inf; end
                end
    
                if ~isfield(v,'h') || ~isfield(v,'T_inf')
                    error('bcs.%s (Convective) must provide h and T_inf (either in .value or as fields).', name);
                end
    
                params.BC.(name).value.h     = hmean;
                params.BC.(name).value.T_inf = tfluid;
    
            case "mixed"
                % Mixed can be problem-specific; safest is to copy bc.value as-is.
                % If you have required subfields (e.g., h, T_inf, q, etc.), validate them here.
                heatPattern(heatPattern ~= 1) = 0;
                params.BC.z1.pttrn = heatPattern;
                heatPattern(heatPattern ~= 1) = 2;
            
                convPattern = heatPattern - 1;
                heatPattern(heatPattern ~= 1) = 0;
                convPattern(2 : end - 1, 2 : end - 1, 2 : end) = 0;
                % Convective parameters (T_inf in Celsius)
                h_matrix = convPattern .* hmean;
                T_inf_matrix = ones(params.GRID.Ny, params.GRID.Nx) .* tfluid; % Ambient air at 22°C
                params.BC.z1.convective_params.h = h_matrix;
                params.BC.z1.convective_params.T_inf = T_inf_matrix;
                
                % Heat flux parameters (q in W/m^2)
                q_applied = qin * 1e0; % 50 kW/m^2
                q_matrix = heatPattern(:, :, 1) .* q_applied;
                params.BC.z1.neumann_params.q = q_matrix;
                
            otherwise
                error('Unknown BC type "%s" in bcs.%s. Allowed: Dirichlet/Neumann/Convective/Mixed.', bc.type, name);
        end
    end

    
    %% 3. Run the Solver
    % T_final_C = solveHeat3D_Euler(params);
    
    % %% 4. Visualize the Results
    % fprintf('Plotting results...\n');
    % 
    % % Create coordinate vectors in mm for plotting
    % y_vec = linspace(0, params.DIMS.Ly, params.GRID.Ny);
    % x_vec = linspace(0, params.DIMS.Lx, params.GRID.Nx);
    % z_vec = linspace(0, params.DIMS.Lz, params.GRID.Nz);
    % 
    % figure('Name', '3D Temperature Field', 'NumberTitle', 'off');
    % [X, Y, Z] = meshgrid(x_vec, y_vec, z_vec);
    % 
    % slice(X, Y, Z, T_final_C, [params.DIMS.Lx/2], [params.DIMS.Ly/2], [params.DIMS.Lz/2]);
    % xlabel('x (mm)');
    % ylabel('y (mm)');
    % zlabel('z (mm)');
    % title(sprintf('Temperature (°C) at t = %.1f s', params.TIME.t_final));
    % cb = colorbar;
    % ylabel(cb, 'Temperature (°C)')
    % axis equal; view(3); grid on;
end