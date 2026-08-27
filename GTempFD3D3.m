function [T, tok] = GTempFD3D3(params, stp)
    tik = tic;
    params.T_initial = permute(params.T_initial, [2 1 3]);
    params.BC.z1.pttrn = permute(params.BC.z1.pttrn, [2 1 3]);
    params.BC.z1.convective_params.h = permute(params.BC.z1.convective_params.h, [2 1 3]);
    params.BC.z1.convective_params.T_inf = permute(params.BC.z1.convective_params.T_inf, [2 1 3]);
    params.BC.z1.neumann_params.q = permute(params.BC.z1.neumann_params.q, [2 1 3]);
    T = heat3d_forwardEuler(params);
    % % FE2DHeatZ1_params  Forward-Euler 2D heat solver for the TOP plane (z1) using your params.mat schema.
    % % Uses in-plane diffusion; pttrn==1 => volumetric generation (q), else convection (h, T_inf).
    % %
    % % INPUT: params (struct) with fields:
    % %   DIMS.Lx, DIMS.Ly, DIMS.Lz
    % %   GRID.Nx, GRID.Ny, GRID.Nz
    % %   MATERIAL.alpha, MATERIAL.k           % alpha = k/(rho*cp)
    % %   TIME.t_final, TIME.dt
    % %   T_initial  (Ny x Nx x Nz)
    % %   BC.z1.type = 'mixed'
    % %   BC.z1.pttrn (Ny x Nx x Nz)           % 1 => generation, 0 => convection
    % %   BC.z1.convective_params.h (Ny x Nx x Nz)
    % %   BC.z1.convective_params.T_inf (Ny x Nx x Nz)
    % %   BC.z1.neumann_params.q (Ny x Nx)     % volumetric heat generation (W/m^3)
    % %
    % % OUTPUT:
    % %   T    : Ny x Nx final temperature at z1 after Nt steps
    % %   meta : diagnostic info (stability, dt_max, etc.)
    % 
    % % ---- pull sizes and steps ----
    % Ny = params.GRID.Ny; Nx = params.GRID.Nx; Nz = params.GRID.Nz;
    % Lx = params.DIMS.Lx; Ly = params.DIMS.Ly; Lz = params.DIMS.Lz;
    % 
    % dx = Lx/(Nx-1);
    % dy = Ly/(Ny-1);
    % dz = Lz/max(Nz,1);            % effective layer thickness for sources/convection
    % 
    % % ---- material/time ----
    % k     = params.MATERIAL.k;
    % alpha = params.MATERIAL.alpha;
    % Cvol  = k / max(alpha, eps);  % rho*cp (per-volume heat capacity)
    % 
    % dt      = params.TIME.dt;
    % t_final = params.TIME.t_final;
    % Nt      = ceil(t_final / dt);
    % 
    % % ---- z1 fields (take top slice along z = Nz) ----
    % iz = 1;                              % top plane
    % Tn = params.T_initial(:,:,iz);
    % pt = squeeze(params.BC.z1.pttrn(:,:,iz));
    % hmap   = squeeze(params.BC.z1.convective_params.h(:,:,iz));
    % T_inf  = squeeze(params.BC.z1.convective_params.T_inf(:,:,iz));
    % q_gen  = params.BC.z1.neumann_params.q;      % already Ny x Nx
    % 
    % % ---- sanity on sizes ----
    % assert(isequal(size(Tn), [Ny, Nx]), 'T_initial top slice must be Ny x Nx.');
    % assert(isequal(size(pt), [Ny, Nx]), 'pttrn top slice must be Ny x Nx.');
    % assert(isequal(size(hmap), [Ny, Nx]), 'h top slice must be Ny x Nx.');
    % assert(isequal(size(T_inf), [Ny, Nx]), 'T_inf top slice must be Ny x Nx.');
    % assert(isequal(size(q_gen), [Ny, Nx]), 'q (generation) must be Ny x Nx.');
    % 
    % isGen  = (pt == 1);
    % isConv = ~isGen;
    % 
    % % ---- stability advisory for explicit 2D FE ----
    % dt_max = (1/(2*alpha)) * 1/(1/dx^2 + 1/dy^2);
    % cfl_ok = (dt <= dt_max + eps);
    % 
    % % ---- precompute coefficients ----
    % % Laplacian term: alpha*dt*∇²T
    % % Source terms:
    % %   Generation:  + dt * q_gen / (rho*cp*dz)
    % %   Convection:  - dt * hmap / (rho*cp*dz) .* (T - T_inf)
    % gen_coef = dt / max(Cvol*dz, eps);
    % conv_coef = dt ./ max(Cvol*dz, eps);   % multiply by h later
    % 
    % % ---- march in time ----
    % for n = 1:Nt
    %     % Mirror edges (Neumann 0-flux): ghost reflection
    %     T_w = cat(2, Tn(:,1),     Tn(:,1:end-1));
    %     T_e = cat(2, Tn(:,2:end), Tn(:,end));
    %     T_n = cat(1, Tn(1,:),     Tn(1:end-1,:));
    %     T_s = cat(1, Tn(2:end,:), Tn(end,:));
    % 
    %     lap = (T_w - 2*Tn + T_e)/dx^2 + (T_n - 2*Tn + T_s)/dy^2;
    % 
    %     % Source field:
    %     S = zeros(Ny, Nx, 'like', Tn);
    %     if any(isGen,'all')
    %         S = S + gen_coef * (q_gen .* double(isGen));
    %     end
    %     if any(isConv,'all')
    %         S = S - conv_coef * (hmap .* (Tn - T_inf) .* double(isConv));
    %     end
    % 
    %     Tn = Tn + alpha*dt*lap + S;
    % end
    % 
    % T = Tn;
    % 
    % % ---- meta ----
    % meta = struct();
    % meta.alpha  = alpha;
    % meta.k      = k;
    % meta.Cvol   = Cvol;
    % meta.dx     = dx;
    % meta.dy     = dy;
    % meta.dz     = dz;
    % meta.dt     = dt;
    % meta.Nt     = Nt;
    % meta.dt_max = dt_max;
    % meta.cfl_ok = logical(cfl_ok);
    % meta.note   = 'Top plane (z1=end). Zero-flux outer edges; generation where pttrn==1 else convection.';
    % 
    % if ~cfl_ok
    %     warning('FE2DHeatZ1_params:CFL', ...
    %         'Forward-Euler may be unstable. dt=%.3g > dt_max≈%.3g', dt, dt_max);
    % end
    % T = cat(3, T, zeros(size(T)));
    % 
    % 
    % SOLVE_HEAT_3D_EULER Solves 3D heat transfer using Forward Euler.
    % 
    %   T = SOLVE_HEAT_3D_EULER(params) takes a structure 'params' loaded from
    %   params.mat and returns the final temperature field T.
    % 
    %   EXPECTED FIELDS in 'params':
    %     - nx, ny, nz: Grid dimensions
    %     - dx, dy, dz: Grid spacing (m)
    %     - dt: Time step (s)
    %     - nt: Number of time steps
    %     - k: Thermal conductivity (W/mK)
    %     - rho: Density (kg/m^3)
    %     - cp: Specific heat (J/kgK)
    %     - T_init: Initial temperature field (3D array)
    %     - pttrn: 2D Mask for z=1 surface (1 = Generation, 0 = Convection)
    %     - q_flux: Heat flux value for generation zones (W/m^2)
    %     - h: Convection coefficient (W/m^2K)
    %     - T_inf: Ambient temperature (K)


%function T_final = heat3d_forwardEuler_from_params(matFile)
%HEAT3D_FORWARDEULER_FROM_PARAMS  3D transient heat conduction (Forward Euler)
%   T_final = heat3d_forwardEuler_from_params("params.mat")
%
% Expects a 1x1 struct named "params" in matFile with (at least):
%   params.DIMS.(Lx,Ly,Lz), params.GRID.(Nx,Ny,Nz)
%   params.MATERIAL.(alpha,k), params.TIME.(dt,t_final)
%   params.T_initial (Nx x Ny x Nz)
%   params.BC with faces: x0,x1,y0,y1,z0 and z1 (mixed)
%
% BC types supported on x0,x1,y0,y1,z0:
%   - "dirichlet": bc.value = T_bc
%   - "neumann"  : bc.value = q  (heat flux into the body) [W/m^2]
%   - "convective": bc.value.h and bc.value.T_inf
%
% z1 is ALWAYS "mixed" with:
%   params.BC.z1.pttrn (Nx x Ny x Nz or Nx x Ny x 1)  uint8 {0,1}
%   params.BC.z1.neumann_params.q (Nx x Ny)           [W/m^2]
%   params.BC.z1.convective_params.h (Nx x Ny x Nz)   [W/m^2-K]
%   params.BC.z1.convective_params.T_inf (Nx x Ny)    [same units as T]
%
% Mixed convention used here:
%   pttrn == 0  -> Neumann (apply heat flux q)
%   pttrn == 1  -> Convective (apply h, T_inf)
% If your convention is inverted, swap the mask condition in makeGhostZ1Mixed().

%     if nargin < 1, matFile = "params.mat"; end
%     S = load(matFile);
%     assert(isfield(S,"params"), 'MAT-file must contain variable "params".');
%     T_final = heat3d_forwardEuler(S.params);
% end

function T = heat3d_forwardEuler(params)

    % --- parse / basic checks ---
    Nx = double(params.GRID.Nx);
    Ny = double(params.GRID.Ny);
    Nz = double(params.GRID.Nz);

    Lx = double(params.DIMS.Lx);
    Ly = double(params.DIMS.Ly);
    Lz = double(params.DIMS.Lz);

    alpha = double(params.MATERIAL.alpha);
    k     = double(params.MATERIAL.k);

    dt      = double(params.TIME.dt);
    t_final = double(params.TIME.t_final);

    T = double(params.T_initial);
    assert(isequal(size(T), [Nx, Ny, Nz]), "T_initial must be Nx x Ny x Nz.");

    dx = Lx/(Nx-1);
    dy = Ly/(Ny-1);
    dz = Lz/(Nz-1);

    % --- stability check (Forward Euler) ---
    dt_crit = 1/(2*alpha*(1/dx^2 + 1/dy^2 + 1/dz^2));
    if dt > dt_crit
        warning("Forward Euler may be unstable: dt=%.3g > dt_crit=%.3g", dt, dt_crit);
    end

    nSteps = max(1, round(t_final/dt));
    BC = params.BC;

    % --- time stepping ---
    for it = 1:nSteps

        % Build 1-cell ghost padding for all 6 faces
        Tpad = zeros(Nx+2, Ny+2, Nz+2);
        Tpad(2:Nx+1, 2:Ny+1, 2:Nz+1) = T;

        % x-faces
        Tpad(1,     2:Ny+1, 2:Nz+1) = makeGhostFace("x0", T, BC.x0, k, dx);
        Tpad(Nx+2,  2:Ny+1, 2:Nz+1) = makeGhostFace("x1", T, BC.x1, k, dx);

        % y-faces
        Tpad(2:Nx+1, 1,     2:Nz+1) = makeGhostFace("y0", T, BC.y0, k, dy);
        Tpad(2:Nx+1, Ny+2,  2:Nz+1) = makeGhostFace("y1", T, BC.y1, k, dy);

        % z0-face
        Tpad(2:Nx+1, 2:Ny+1, 1)     = makeGhostFace("z0", T, BC.z0, k, dz);

        % z1-face (mixed, always)
        Tpad(2:Nx+1, 2:Ny+1, Nz+2)  = makeGhostZ1Mixed(T, BC.z1, k, dz);

        % Laplacian on all nodes (including boundaries, via ghosts)
        Tc  = Tpad(2:Nx+1, 2:Ny+1, 2:Nz+1);
        Txx = (Tpad(3:Nx+2, 2:Ny+1, 2:Nz+1) - 2*Tc + Tpad(1:Nx,   2:Ny+1, 2:Nz+1)) / dx^2;
        Tyy = (Tpad(2:Nx+1, 3:Ny+2, 2:Nz+1) - 2*Tc + Tpad(2:Nx+1, 1:Ny,   2:Nz+1)) / dy^2;
        Tzz = (Tpad(2:Nx+1, 2:Ny+1, 3:Nz+2) - 2*Tc + Tpad(2:Nx+1, 2:Ny+1, 1:Nz  )) / dz^2;

        T = Tc + alpha * dt * (Txx + Tyy + Tzz);

        % Re-enforce Dirichlet faces exactly (if any changed in input files)
        T = enforceDirichlet(T, BC);
    end
end

function G = makeGhostFace(face, T, bc, k, d)
% Returns the ghost plane values (same size as the boundary plane in T)

    bcType = lower(string(bc.type));

    switch face
        case "x0"
            T_in   = T(2,:,:);
            T_bdry = T(1,:,:);
            isMin  = true;
        case "x1"
            T_in   = T(end-1,:,:);
            T_bdry = T(end,:,:);
            isMin  = false;
        case "y0"
            T_in   = T(:,2,:);
            T_bdry = T(:,1,:);
            isMin  = true;
        case "y1"
            T_in   = T(:,end-1,:);
            T_bdry = T(:,end,:);
            isMin  = false;
        case "z0"
            assert(size(T,3) >= 2, "Nz must be >= 2 to use z0 ghosting.");
            T_in   = T(:,:,2);
            T_bdry = T(:,:,1);
            isMin  = true;
        otherwise
            error("Unsupported face: %s", face);
    end

    % Ensure T_in, T_bdry are double arrays
    T_in   = double(T_in);
    T_bdry = double(T_bdry);

    if bcType == "dirichlet"
        Tbc = double(bc.value);
        G = 2*Tbc - T_in;  % second-order consistent ghost

    elseif bcType == "neumann"
        q = getNeumannFlux(bc); % [W/m^2] into the body (positive adds heat)
        q = expandToSize(q, size(T_in));

        % Convert flux to gradient along +axis: dT/dx = sgn*q/k  (min:+, max:-)
        sgn = +1;
        if ~isMin, sgn = -1; end
        g = sgn * (q ./ k); % [K/m]

        if isMin
            % (T_in - Tghost)/(2d) = g  -> Tghost = T_in - 2d*g
            G = T_in - 2 * d .* g;
        else
            % (Tghost - T_in)/(2d) = g -> Tghost = T_in + 2d*g
            G = T_in + 2 * d .* g;
        end

    elseif bcType == "convective"
        [h, Tinf] = getConvectiveParams(bc);
        h    = expandToSize(h,    size(T_in));
        Tinf = expandToSize(Tinf, size(T_in));

        % Use centered derivative at the boundary:
        % min face:  k*dT/dx = h(Ts-Tinf)  -> (T_in - Tghost)/(2d) = (h/k)(Ts-Tinf)
        % max face: -k*dT/dx = h(Ts-Tinf)  -> (Tghost - T_in)/(2d) = -(h/k)(Ts-Tinf)
        % Both yield: Tghost = T_in - 2d*(h/k)*(Ts - Tinf)
        G = T_in - 2*d .* (h./k) .* (T_bdry - Tinf);

    else
        error("Unsupported BC type '%s' on face '%s'.", bcType, face);
    end
end

function G = makeGhostZ1Mixed(T, z1bc, k, dz)
% Mixed BC at z1 (top face, k=Nz). Returns ghost plane for z = Lz + dz.

    assert(lower(string(z1bc.type)) == "mixed", "BC.z1.type must be 'mixed'.");

    Nx = size(T,1); Ny = size(T,2); Nz = size(T,3);
    assert(Nz >= 2, "Nz must be >= 2 for z1 mixed BC (needs interior neighbor).");

    % Use the top slice in the pttrn/h arrays
    pt = z1bc.pttrn;
    if ndims(pt) == 3
        maskConv = (pt(:, :, 1) == 0);  % 1 -> convective, 0 -> neumann (default here)
    else
        maskConv = (pt(:, :) == 0);
    end

    q    = double(z1bc.neumann_params.q);              % Nx x Ny
    h3   = double(z1bc.convective_params.h);           % Nx x Ny x Nz
    h    = h3(:, :, 1);                                 % Nx x Ny
    Tinf = double(z1bc.convective_params.T_inf);       % Nx x Ny

    T_in   = double(T(:, :, Nz - 1));  % neighbor inside
    T_surf = double(T(:, :, Nz));    % boundary plane

    % Neumann part (z1 is a max face): dT/dz = -q/k
    ghostNeu = T_in + 2 * dz .* (q ./ k);

    % Convective part: ghost = T_in - 2dz*(h/k)*(Ts - Tinf)
    ghostCon = T_in - 2 * dz .* (h ./ k) .* (T_surf - Tinf);

    G2 = ghostNeu;
    G2(maskConv) = ghostCon(maskConv);

    G = reshape(G2, [Nx, Ny, 1]); % pad assignment expects Nx x Ny x 1
end

function T = enforceDirichlet(T, BC)
% Clamp any faces that are set to Dirichlet

    if isfield(BC,"x0") && lower(string(BC.x0.type)) == "dirichlet"
        T(1,:,:) = double(BC.x0.value);
    end
    if isfield(BC,"x1") && lower(string(BC.x1.type)) == "dirichlet"
        T(end,:,:) = double(BC.x1.value);
    end
    if isfield(BC,"y0") && lower(string(BC.y0.type)) == "dirichlet"
        T(:,1,:) = double(BC.y0.value);
    end
    if isfield(BC,"y1") && lower(string(BC.y1.type)) == "dirichlet"
        T(:,end,:) = double(BC.y1.value);
    end
    if isfield(BC,"z0") && lower(string(BC.z0.type)) == "dirichlet"
        T(:,:,1) = double(BC.z0.value);
    end
    if isfield(BC,"z1") && isfield(BC.z1,"type") && lower(string(BC.z1.type)) == "dirichlet"
        T(:,:,end) = double(BC.z1.value);
    end
end

function q = getNeumannFlux(bc)
% Returns q (heat flux into the body) [W/m^2]
    if isstruct(bc.value)
        if isfield(bc.value,"q")
            q = double(bc.value.q);
        else
            error("Neumann BC.value is a struct but has no field 'q'.");
        end
    else
        q = double(bc.value);
    end
end

function [h, Tinf] = getConvectiveParams(bc)
    if isstruct(bc.value)
        assert(isfield(bc.value,"h") && isfield(bc.value,"T_inf"), ...
            "Convective BC.value must contain fields 'h' and 'T_inf'.");
        h    = double(bc.value.h);
        Tinf = double(bc.value.T_inf);
    else
        error("Convective BC.value must be a struct with fields 'h' and 'T_inf'.");
    end
end

function A = expandToSize(A, targetSize)
% Expands scalar / 2D to targetSize by replication (no interpolation)
    if isscalar(A)
        A = repmat(A, targetSize);
        return;
    end
    sz = size(A);
    if isequal(sz, targetSize)
        return;
    end
    % allow 2D -> 3D when third dim is singleton
    if numel(sz) == 2 && numel(targetSize) == 3 && targetSize(3) == 1
        if isequal(sz, targetSize(1:2))
            A = reshape(A, [targetSize(1), targetSize(2), 1]);
            return;
        end
    end
    error("Cannot expand array of size [%s] to [%s].", num2str(sz), num2str(targetSize));
end

T = permute(T, [2 1 3]);
tok = sprintf('Elapsed: %.6f s', toc(tik));

end
