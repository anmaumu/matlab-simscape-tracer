function r = tracer_net(P, varargin)
%TRACER_NET MATLAB prototype of the extended mylib tracer network (same equations as the .ssc files).
%   Chain of volumes V(1..N) linked by elements elem{1..N+1}:
%     inlet boundary -elem{1}- V1 -elem{2}- V2 ... VN -elem{N+1}- outlet boundary
%   Element types:
%     'orifice' : Local Restriction (IL) equations (Re_cr based dp_cr, area ratio, pressure recovery)
%     'pipe'    : Pipe friction (laminar / smoothed transition / turbulent Haaland), mdot is algebraic
%     'source'  : fixed mass flow from the inlet boundary (only as elem{1})
%     'closed'  : no flow
%   Density depends on pressure: rho = rho0*exp((p - pRef)/beta).
%   Tracer flow: upwind advection + diffusion rho*Dc*S/L*(cA - cB).
%   Solved as an index-1 DAE with ode15s (pipe mass flows are algebraic states).
%
%   Element equations can be called directly for testing, e.g.
%   m  = tracer_net('orificeFlow', dp, rho, E, F);
%   dp = tracer_net('pipeDp', m, rho, E, F);

    if ischar(P)
        fn = str2func(P);
        r = fn(varargin{:});
        return
    end

    N  = numel(P.V);
    F  = P.fluid;
    types = cellfun(@(e) e.type, P.elem, 'UniformOutput', false);
    iPipe = find(strcmp(types, 'pipe'));
    nP = numel(iPipe);

    x0 = initialState(P, N, iPipe);
    M = blkdiag(eye(2*N), zeros(nP));
    opts = odeset('Mass', M, 'MassSingular', 'yes', 'RelTol', 1e-7, ...
        'AbsTol', [1e-2*ones(N, 1); 1e-10*ones(N, 1); 1e-10*ones(nP, 1)]);

    edges = unique([0, P.tOn, P.tOff, P.tEnd]);
    t = []; x = [];
    for k = 1:numel(edges) - 1
        [tk, xk] = ode15s(@(tt, xx) rhs(tt, xx, P, N), [edges(k), edges(k+1)], x0, opts);
        t = [t; tk]; %#ok<AGROW>
        x = [x; xk]; %#ok<AGROW>
        x0 = xk(end, :).';
    end

    nT = numel(t);
    m = zeros(nT, N + 1);
    mc = zeros(nT, N + 1);
    for k = 1:nT
        [m(k, :), mc(k, :)] = network(x(k, :).', P, N);
    end
    p = x(:, 1:N);
    c = x(:, N+1:2*N);
    rho = density(p, F);
    inj = P.mInj * (t >= P.tOn & t < P.tOff);

    r.t = t; r.p = p; r.c = c; r.rho = rho; r.mdot = m; r.mdot_c = mc; r.params = P;
    r.massInSystem   = rho * P.V(:);
    r.tracerInSystem = (rho .* c) * P.V(:);
    r.massBalanceErr   = (r.massInSystem(end) - r.massInSystem(1)) ...
        - trapz(t, m(:, 1) - m(:, end) + inj);
    r.tracerBalanceErr = (r.tracerInSystem(end) - r.tracerInSystem(1)) ...
        - trapz(t, mc(:, 1) - mc(:, end) + inj * P.cInj);
    r.injectedTracer = P.mInj * P.cInj * (P.tOff - P.tOn);
end

function dx = rhs(t, x, P, N)
    F = P.fluid;
    p = x(1:N);
    c = x(N+1:2*N);
    [m, mc, res] = network(x, P, N);
    rho = density(p, F);

    inj = zeros(N, 1);
    if t >= P.tOn && t < P.tOff
        inj(P.injVol) = P.mInj;
    end
    netM  = m(1:N)  - m(2:N+1)  + inj;
    netMc = mc(1:N) - mc(2:N+1) + inj * P.cInj;

    dp = F.beta ./ (rho .* P.V(:)) .* netM;            % V*(drho/dp)*p' = mdot, drho/dp = rho/beta
    dc = (netMc - c .* netM) ./ (rho .* P.V(:));       % rho*V*c' = mdot_c - c*mdot
    dx = [dp; dc; res];
end

function [m, mc, res] = network(x, P, N)
    F = P.fluid;
    pAll = [P.pIn; x(1:N); P.pOut];
    cAll = [P.cIn; x(N+1:2*N); P.cOut];
    mPipe = x(2*N+1:end);
    rhoAll = density(pAll, F);

    m = zeros(N + 1, 1);
    mc = zeros(N + 1, 1);
    res = zeros(numel(mPipe), 1);
    j = 0;
    for k = 1:N + 1
        E = P.elem{k};
        a = k; b = k + 1;
        rhoAvg = (rhoAll(a) + rhoAll(b)) / 2;
        dp = pAll(a) - pAll(b);
        switch E.type
            case 'source'
                m(k) = E.mdot;
                mc(k) = E.mdot * cAll(a);
                continue
            case 'closed'
                continue
            case 'orifice'
                m(k) = orificeFlow(dp, rhoAvg, E, F);
                Sd = E.SR;
                Ld = sqrt(4 * E.SR / pi);
            case 'pipe'
                j = j + 1;
                m(k) = mPipe(j);
                res(j) = dp - pipeDp(m(k), rhoAvg, E, F);
                Sd = pi * E.D^2 / 4;
                Ld = E.L;
        end
        w = (1 + tanh(m(k) / P.mEps)) / 2;
        cUp = w * cAll(a) + (1 - w) * cAll(b);
        G = rhoAvg * F.Dc * Sd / Ld;                    % Diffusive conductance [kg/s]
        mc(k) = m(k) * cUp + G * (cAll(a) - cAll(b));
    end
end

function x0 = initialState(P, N, iPipe)
    F = P.fluid;
    E1 = P.elem{1};
    EN = P.elem{end};
    if strcmp(EN.type, 'closed')
        p0 = P.pOut * ones(N, 1);
    elseif strcmp(E1.type, 'source')
        % March upstream from the outlet with the source flow (steady state)
        p0 = zeros(N, 1);
        pDown = P.pOut;
        for k = N + 1:-1:2
            p0(k-1) = pDown + elementDp(E1.mdot, P.elem{k}, density(pDown, F), F);
            pDown = p0(k-1);
        end
    else
        pLin = linspace(P.pIn, P.pOut, N + 2).';
        p0 = pLin(2:end-1);
    end
    if isfield(P, 'c0')
        c0 = P.c0(:);
    else
        c0 = zeros(N, 1);
    end

    % Consistent algebraic pipe flows for the initial pressures
    pAll = [P.pIn; p0; P.pOut];
    mP = zeros(numel(iPipe), 1);
    for j = 1:numel(iPipe)
        k = iPipe(j);
        dp = pAll(k) - pAll(k+1);
        rhoAvg = mean(density(pAll(k:k+1), F));
        mP(j) = fzero(@(mm) pipeDp(mm, rhoAvg, P.elem{k}, F) - dp, 0);
    end
    x0 = [p0; c0; mP];
end

function dp = elementDp(mdot, E, rho, F)
    switch E.type
        case 'pipe'
            dp = pipeDp(mdot, rho, E, F);
        case 'orifice'
            dp = fzero(@(d) orificeFlow(d, rho, E, F) - mdot, [-1e9, 1e9]);
        otherwise
            dp = 0;
    end
end

function rho = density(p, F)
    rho = F.rho0 * exp((p - F.pRef) / F.beta);
end

function m = orificeFlow(dp, rho, E, F)
    % Local Restriction (IL)
    dpCr  = pi/4 * rho / (2 * E.SR) * (E.Recr * F.nu / E.Cd)^2;
    ratio = E.SR / E.S;
    root  = sqrt(1 - ratio^2 * (1 - E.Cd^2));
    PRrec = (root - E.Cd * ratio) / (root + E.Cd * ratio);
    PR    = E.recovery * PRrec + (1 - E.recovery);
    m = E.Cd * E.SR * sqrt(2 * rho / (PR * (1 - ratio^2))) * dp / (dp^2 + dpCr^2)^0.25;
end

function dp = pipeDp(m, rho, E, F)
    % Pipe friction: laminar (Hagen-Poiseuille) / smoothed transition / turbulent (Haaland)
    S  = pi * E.D^2 / 4;
    Re = abs(m) * E.D / (S * rho * F.nu);
    dpLam = E.lam * F.nu * (E.L + E.Leq) * m / (2 * E.D^2 * S);
    f = haaland(max(Re, E.ReL), E.rough / E.D);
    dpTur = f * (E.L + E.Leq) * m * abs(m) / (2 * rho * E.D * S^2);
    s = smoothstep((Re - E.ReL) / (E.ReT - E.ReL));
    dp = (1 - s) * dpLam + s * dpTur;
end

function f = haaland(Re, relRough)
    f = 1 / (-1.8 * log10(6.9 / Re + (relRough / 3.7)^1.11))^2;
end

function s = smoothstep(z)
    z = min(max(z, 0), 1);
    s = 3*z^2 - 2*z^3;
end
