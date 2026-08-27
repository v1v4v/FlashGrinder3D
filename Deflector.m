function [state_out, history, time] = Deflector(state_in, F_tip, support, dt_global, params)
tik = tic;
% BEAM_TIP_VIBRATION_STEP (2-DOF: x & z)
% Linear cantilever beam with fixed base. The TIP is attached to a spring-damper
% to ground. Dynamics are in two lateral directions, x and z, both ⟂ to shaft axis.
% The beam contributes its own stiffness/damping in each direction. The external
% tip force is constant during this call (substepping), but can change next call.
% Time integration: Newmark-β (β=1/4, γ=1/2), unconditionally stable.
%
% Inputs
%   state_in.x, state_in.v, (optional) state_in.a : 2x1 vectors [x; z] at start
%   F_tip           : 2x1 vector [Fx; Fz] (constant during this call)
%   support.k       : scalar | [kx kz] | 2x2 matrix (support stiffness)
%   support.c       : scalar | [cx cz] | 2x2 matrix (support damping)
%   dt_global       : global time step (this call spans this duration)
%   params : struct (key fields below; others optional)
%       Mass:
%         .m        : scalar | [mx mz] | 2x2 (effective tip mass matrix)
%         OR provide .rho,.A,.L (and optional .alpha_m ~ 0.236) to estimate m.
%       Beam stiffness:
%         .k_beam   : scalar | [kx kz] | 2x2
%         OR provide .E,.L and either { .Ix,.Iz } or single .I (Ix=Iz=I) to set
%         kx = 3*E*Iz/L^3 (deflection in x bends about z), kz = 3*E*Ix/L^3.
%       Beam damping:
%         .c_beam   : scalar | [cx cz] | 2x2. If omitted and .zeta_beam is
%         scalar|[ζx ζz] AND M,K are diagonal, sets c = 2 ζ ⊙ sqrt(k ⊙ m).
%       Integration / substepping:
%         .n_substeps (default 10), .beta (1/4), .gamma (1/2)
%
% Outputs
%   state_out.x, state_out.v, state_out.a : 2x1 tip state at end of this call
%   state_out.K, state_out.C, state_out.M : matrices used this call
%   history (optional): .t_local (1×(n+1)), .x,.v,.a ((n+1)×2)
%
% Example:
%   params.E=210e9; params.L=0.25; params.Ix=1.2e-8; params.Iz=8.0e-9;
%   params.rho=7800; params.A=2.0e-4; params.n_substeps=20; params.zeta_beam=[0.01 0.01];
%   state = struct('x',[0;0],'v',[0;0]);
%   [state,~] = beam_tip_vibration_step(state,[5;12], 2e5, 150, 1e-3, params);

    % ---- defaults ----
    beta  = getOr(params,'beta',  1 / 4);
    gamma = getOr(params,'gamma', 1 / 2);
    n     = getOr(params,'n_substeps', 10);
    dt    = dt_global / n;

    % ---- beam stiffness K_beam ----
    if ~isfield(params,'k_beam') || isempty(params.k_beam)
        need = all(isfield(params,{'E','L'}));
        if need && isfield(params,'Ix') && isfield(params,'Iz')
            kx = params.E * params.A / params.L; % deflection in x ⇒ bend about z (Iz)
            kz = params.E * params.A / params.L; % deflection in z ⇒ bend about x (Ix)
            Kb = diag([kx kz]);
        elseif need && isfield(params,'I')
            k = params.E * params.A / params.L;
            Kb = k * eye(2);
        else
            error('Provide params.k_beam or (E,L and Ix/Iz or I).');
        end
    else
        Kb = to2x2(params.k_beam,'params.k_beam');
    end

    % ---- mass M ----
    if ~isfield(params,'m') || isempty(params.m)
        if all(isfield(params,{'rho','A','L'}))
            alpha_m = getOr(params,'alpha_m',0.236);
            m = alpha_m * params.rho * params.A * params.L; % per lateral DOF
            M = m*eye(2);
        else
            error('Provide params.m or (rho,A,L).');
        end
    else
        M = to2x2(params.m,'params.m');
    end

    % ---- beam damping C_beam ----
    if ~isfield(params,'c_beam') || isempty(params.c_beam)
        if isfield(params,'zeta_beam') && ~isempty(params.zeta_beam) && isdiaglike(M) && isdiaglike(Kb)
            zeta = to2vec(params.zeta_beam,'params.zeta_beam');
            mdiag = diag(M);
            kdiag = diag(Kb);
            cdiag = 2 .* zeta .* sqrt(kdiag .* mdiag);
            Cb = diag(cdiag);
        else
            Cb = zeros(2);
        end
    else
        Cb = to2x2(params.c_beam,'params.c_beam');
    end

    % ---- support K,C (per-call) ----
    % Ks = to2x2(support.k,'support.k');
    % Cs = to2x2(support.c,'support.c');

    % Totals this call
    K = Kb .* support;
    C = Cb .* support;

    % ---- initial state ----
    x = requireVec2(state_in,'x');  % 2x1
    v = requireVec2(state_in,'v');  % 2x1
    if isfield(state_in,'a') && ~isempty(state_in.a)
        a = requireVec2(state_in,'a');
    else
        F = to2vec(F_tip,'F_tip');
        a = M \ (F - C * v - K * x);
    end

    % ---- history (optional) ----
    if nargout > 1
        history.t_local = (0:n)*dt;
        history.x = zeros(n+1,2);
        history.v = history.x;
        history.a = history.x;

        history.x(1,:) = x.';
        history.v(1,:) = v.';
        history.a(1,:) = a.';
    end

    % ---- Newmark-β substepping ----
    F = to2vec(F_tip,'F_tip');  % constant within this call
    Mhat = []; % lazily formed each step in case K/C change (they don't here, but cheap)
    for i = 1:n
        x_pred = x + dt*v + dt^2*(0.5 - beta)*a;
        v_pred = v + dt*(1 - gamma)*a;

        if isempty(Mhat)
            Mhat = M + gamma*dt*C + beta*dt^2*K;
        end
        rhs  = F - C*v_pred - K*x_pred;

        a_new = Mhat \ rhs;
        x_new = x_pred + beta*dt^2*a_new;
        v_new = v_pred + gamma*dt*a_new;

        x = x_new; v = v_new; a = a_new;

        if nargout > 1
            history.x(i+1,:) = x.'; history.v(i+1,:) = v.'; history.a(i+1,:) = a.';
        end
    end

    % ---- outputs ----
    state_out.x = x;
    state_out.v = v;
    state_out.a = a;

    tok = toc(tik);
    time = ['The elapsed time for computing the tool deflection is ', num2str(tok), ' seconds'];

% main function
end

    % ----------------- helpers -----------------
    function val = getOr(s,fld,defaultVal)
        if isfield(s,fld) && ~isempty(s.(fld))
            val = s.(fld);
        else
            val = defaultVal;
        end
    end
    
    function X = to2x2(A,name)
        if isscalar(A)
            X = A*eye(2);
        elseif isvector(A) && numel(A)==2
            A = A(:);
            X = diag(A);
        elseif isequal(size(A),[2 2])
            X = A;
        else
            error('%s must be scalar, 1x2/2x1, or 2x2.', name);
        end
    end
    
    function v = to2vec(a,name)
        if isvector(a) && numel(a)==2
            v = a(:);
        elseif isscalar(a)
            error('%s must be a 2x1 vector [ax; az] for 2-DOF mode.', name);
        else
            error('%s must be scalar or 2-vector.', name);
        end
    end
    
    function v = requireVec2(s, fld)
        if ~isfield(s,fld) || isempty(s.(fld))
            v = [0;0];
        else
            v = to2vec(s.(fld), sprintf('state_in.%s',fld));
        end
    end
    
    function tf = isdiaglike(M)
        tf = isequal(size(M),[2 2]) && abs(M(1,2))<eps && abs(M(2,1))<eps;
    end