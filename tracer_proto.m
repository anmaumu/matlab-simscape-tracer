function results = tracer_proto(varargin)
%TRACER_PROTO MATLAB prototype of the mylib tracer network (same equations as the .ssc files).
%   Inlet - [Volume - Orifice] x N - Reservoir(c=0)
%   inletMode 'pressure': inlet is a reservoir + orifice (fixed pressure)
%   inletMode 'flow'    : inlet is a fixed mass flow source (pump)
%   Pure tracer (c=1) is injected into volume 1 during [tOn, tOff].
%
%   Any parameter can be overridden by name-value pairs, e.g.
%   r = tracer_proto('N', 30, 'inletMode', 'flow', 'showPlot', false);

    % Fluid and geometry (SI units)
    P.rho  = 850;        % Density [kg/m^3]
    P.beta = 1.5e9;      % Bulk modulus [Pa]
    P.Cd   = 0.64;       % Discharge coefficient
    P.area = 1e-5;       % Orifice area [m^2]
    P.dpCr = 1e3;        % Transition pressure difference [Pa]
    P.mEps = 1e-6;       % Upwind smoothing width [kg/s]
    P.N    = 10;         % Number of volume segments
    P.Vtot = 1e-3;       % Total volume [m^3]

    % Boundary conditions
    P.inletMode = 'pressure';
    P.pIn  = 0.3e6;  P.cIn  = 0;   % Inlet reservoir (pressure mode)
    P.mInFix = 0.03555;            % Inlet mass flow (flow mode) [kg/s]
    P.pOut = 0.1e6;  P.cOut = 0;   % Outlet reservoir
    P.mInj = 0.005;  P.cInj = 1;   % Injector [kg/s], concentration
    P.tOn  = 10;     P.tOff = 70;
    P.tEnd = 150;
    P.showPlot = true;
    P.verbose  = true;

    for k = 1:2:numel(varargin)
        P.(varargin{k}) = varargin{k+1};
    end
    P.V = P.Vtot / P.N;
    N = P.N;

    x0 = [initialPressure(P); zeros(N, 1)];

    % Integrate piecewise so the injector on/off edges are not smeared
    opts = odeset('RelTol', 1e-7, 'AbsTol', [1e-2*ones(N, 1); 1e-10*ones(N, 1)]);
    edges = [0, P.tOn, P.tOff, P.tEnd];
    t = []; x = [];
    for k = 1:numel(edges) - 1
        [tk, xk] = ode15s(@(tt, xx) rhs(tt, xx, P), [edges(k), edges(k+1)], x0, opts);
        t = [t; tk]; %#ok<AGROW>
        x = [x; xk]; %#ok<AGROW>
        x0 = xk(end, :).';
    end

    % Post-processing
    nT = numel(t);
    mFlow = zeros(nT, N + 1);
    mcFlow = zeros(nT, N + 1);
    for k = 1:nT
        [mFlow(k, :), mcFlow(k, :)] = flows(x(k, :).', P);
    end
    p = x(:, 1:N);
    c = x(:, N+1:end);

    injected = P.mInj * P.cInj * (P.tOff - P.tOn);
    exited   = trapz(t, mcFlow(:, end)) - trapz(t, min(mcFlow(:, 1), 0));
    stored   = P.rho * P.V * sum(c(end, :));

    if P.verbose
        [~, iSS] = min(abs(t - P.tOff));
        cTheory = P.mInj / mFlow(iSS, end);
        fprintf('Flow before injection : in %.5f / out %.5f kg/s\n', mFlow(find(t < P.tOn, 1, 'last'), [1 end]));
        fprintf('Flow during injection : in %.5f / out %.5f kg/s\n', mFlow(iSS, [1 end]));
        fprintf('Outlet c at t=%g s    : %.5f (theory %.5f)\n', P.tOff, c(iSS, end), cTheory);
        fprintf('Tracer injected %.5f kg, exited %.5f kg, stored %.5f kg, error %.2e kg\n', ...
            injected, exited, stored, injected - exited - stored);
    end

    if P.showPlot
        cTheory = P.mInj / mFlow(find(t <= P.tOff, 1, 'last'), end);
        fig = figure('Name', 'Tracer prototype', 'Position', [100 100 800 700]);
        subplot(3, 1, 1);
        plot(t, mFlow(:, 1), t, mFlow(:, end), 'LineWidth', 1.2);
        xline([P.tOn P.tOff], ':');
        ylabel('mdot [kg/s]'); legend('inlet', 'outlet', 'Location', 'best'); grid on;
        title('Total mass flow');
        subplot(3, 1, 2);
        plot(t, c(:, [1 round(N/2) N]), 'LineWidth', 1.2);
        xline([P.tOn P.tOff], ':'); yline(cTheory, '--');
        ylabel('c [-]'); legend('vol 1', sprintf('vol %d', round(N/2)), sprintf('vol %d', N), 'Location', 'best'); grid on;
        title('Tracer mass fraction');
        subplot(3, 1, 3);
        plot(t, p(:, [1 N]) / 1e6, 'LineWidth', 1.2);
        xline([P.tOn P.tOff], ':');
        ylabel('p [MPa]'); xlabel('t [s]'); legend('vol 1', sprintf('vol %d', N), 'Location', 'best'); grid on;
        title('Pressure');
        exportgraphics(fig, fullfile(fileparts(mfilename('fullpath')), 'tracer_proto.png'));
    end

    results = struct('t', t, 'p', p, 'c', c, 'mdot', mFlow, 'mdot_c', mcFlow, ...
        'injected', injected, 'exited', exited, 'stored', stored, 'params', P);
end

function p0 = initialPressure(P)
    % Steady pressure profile without injection
    N = P.N;
    if strcmp(P.inletMode, 'flow')
        k = P.Cd * P.area * sqrt(2 * P.rho);
        dpSS = fzero(@(dp) k * dp / (dp^2 + P.dpCr^2)^0.25 - P.mInFix, [0, 1e8]);
        p0 = P.pOut + (N:-1:1).' * dpSS;
    else
        pAll = linspace(P.pIn, P.pOut, N + 2).';
        p0 = pAll(2:end-1);
    end
end

function dx = rhs(t, x, P)
    N = P.N;
    c = x(N+1:end);
    [m, mc] = flows(x, P);

    inj = zeros(N, 1);
    if t >= P.tOn && t < P.tOff
        inj(1) = P.mInj;
    end

    netM  = m(1:N).'  - m(2:N+1).'  + inj;          % Mass flow into each volume
    netMc = mc(1:N).' - mc(2:N+1).' + inj * P.cInj; % Tracer flow into each volume

    dp = P.beta / (P.rho * P.V) * netM;             % V*rho/beta * p' = mdot
    dc = (netMc - c .* netM) / (P.rho * P.V);       % rho*V * c' = mdot_c - c*mdot
    dx = [dp; dc];
end

function [m, mc] = flows(x, P)
    % Flow k enters volume k (k = 1: inlet), flow N+1 leaves to the outlet reservoir
    N = P.N;
    pAll = [P.pIn; x(1:N); P.pOut];
    cAll = [P.cIn; x(N+1:end); P.cOut];
    dp = pAll(1:end-1) - pAll(2:end);
    m = P.Cd * P.area * sqrt(2 * P.rho) * dp ./ (dp.^2 + P.dpCr^2).^0.25;
    if strcmp(P.inletMode, 'flow')
        m(1) = P.mInFix;
    end
    w = (1 + tanh(m / P.mEps)) / 2;
    mc = m .* (w .* cAll(1:end-1) + (1 - w) .* cAll(2:end));
    m = m.';
    mc = mc.';
end
