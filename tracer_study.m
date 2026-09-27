%% TRACER_STUDY  Parameter studies with tracer_proto
outDir = fileparts(mfilename('fullpath'));
common = {'showPlot', false, 'verbose', false};

%% 1. Inlet boundary: fixed pressure vs fixed flow (step injection 0.005 kg/s)
modes = {'pressure', 'flow'};
fig1 = figure('Name', 'Inlet boundary', 'Position', [100 100 900 600]);
fprintf('--- 1. Inlet boundary (step injection %.3f kg/s) ---\n', 0.005);
fprintf('%-9s %12s %12s %12s %12s %10s\n', 'inlet', 'in before', 'in during', 'out during', 'out incr', 'c_out');
for k = 1:numel(modes)
    r = tracer_proto('inletMode', modes{k}, common{:});
    iBefore = find(r.t < r.params.tOn, 1, 'last');
    iDuring = find(r.t <= r.params.tOff, 1, 'last');
    fprintf('%-9s %12.5f %12.5f %12.5f %12.5f %10.5f\n', modes{k}, r.mdot(iBefore, 1), ...
        r.mdot(iDuring, 1), r.mdot(iDuring, end), r.mdot(iDuring, end) - r.mdot(iBefore, end), r.c(iDuring, end));
    subplot(2, 2, k);
    plot(r.t, r.mdot(:, 1), r.t, r.mdot(:, end), 'LineWidth', 1.2); grid on;
    title(['mdot, inlet = ' modes{k}]); ylabel('kg/s'); legend('inlet', 'outlet', 'Location', 'best');
    subplot(2, 2, k + 2);
    plot(r.t, r.c(:, end), 'LineWidth', 1.2); grid on;
    title(['outlet c, inlet = ' modes{k}]); xlabel('t [s]');
end
exportgraphics(fig1, fullfile(outDir, 'study1_inlet.png'));

%% 2. Segment count N (flow inlet, step injection)
Ns = [1 3 10 30 100];
fig2 = figure('Name', 'Segment count', 'Position', [100 100 800 450]);
hold on; grid on;
fprintf('\n--- 2. Segment count (flow inlet, step injection at t = 10 s) ---\n');
fprintf('%5s %14s %14s\n', 'N', 't(c=10%)-10', 't(c=90%)-10');
for N = Ns
    r = tracer_proto('inletMode', 'flow', 'N', N, 'tOff', 149, 'tEnd', 150, common{:});
    cFinal = r.params.mInj / (r.params.mInFix + r.params.mInj);
    cOut = r.c(:, end) / cFinal;
    % Interpolate on the rising part (solver output can be sparse)
    iRise = find(r.t > r.params.tOn & r.t < r.params.tOff);
    [cU, iU] = unique(cOut(iRise));
    tRise = r.t(iRise(iU));
    t10 = interp1(cU, tRise, 0.1) - r.params.tOn;
    t90 = interp1(cU, tRise, 0.9) - r.params.tOn;
    fprintf('%5d %14.2f %14.2f\n', N, t10, t90);
    plot(r.t - r.params.tOn, cOut, 'LineWidth', 1.2, 'DisplayName', sprintf('N = %d', N));
end
tauInj = r.params.rho * r.params.Vtot / (r.params.mInFix + r.params.mInj);
xline(tauInj, '--', 'DisplayName', sprintf('\\tau = %.1f s', tauInj));
xlim([0 100]); xlabel('t - t_{on} [s]'); ylabel('c_{out} / c_{final}');
title('Outlet response to step injection'); legend('Location', 'southeast');
exportgraphics(fig2, fullfile(outDir, 'study2_segments.png'));

%% 3. Pulse injection -> residence time distribution (RTD)
tPulse = 1;
tau = r.params.rho * r.params.Vtot / r.params.mInFix;
fig3 = figure('Name', 'RTD', 'Position', [100 100 800 450]);
hold on; grid on;
fprintf('\n--- 3. RTD (flow inlet, %.0f s pulse), tau = %.2f s ---\n', tPulse, tau);
fprintf('%5s %10s %10s %14s %14s %12s\n', 'N', 'mean', 'tau', 'var/tau^2', 'theory 1/N', 'recovered');
colors = lines(numel(Ns));
for k = 1:numel(Ns)
    N = Ns(k);
    r = tracer_proto('inletMode', 'flow', 'N', N, 'tOn', 10, 'tOff', 10 + tPulse, ...
        'tEnd', 400, common{:});
    ts = r.t - (r.params.tOn + tPulse/2);            % time from pulse centre
    E  = r.mdot_c(:, end) / r.injected;              % RTD E(t) [1/s]
    tm = trapz(ts, ts .* E);
    s2 = trapz(ts, (ts - tm).^2 .* E);
    fprintf('%5d %10.2f %10.2f %14.4f %14.4f %12.4f\n', N, tm, tau, s2 / tau^2, 1 / N, trapz(ts, E));
    plot(ts, E, 'LineWidth', 1.2, 'Color', colors(k, :), 'DisplayName', sprintf('N = %d', N));
    % Tanks-in-series theory
    tt = linspace(0, 100, 500);
    ti = tau / N;
    Eth = exp((N-1) * log(tt / ti) - tt / ti - gammaln(N)) / ti;
    plot(tt, Eth, '--', 'Color', colors(k, :), 'HandleVisibility', 'off');
end
xlim([0 100]); xlabel('t [s]'); ylabel('E(t) [1/s]');
title('RTD: simulation (solid) vs tanks-in-series theory (dashed)'); legend;
exportgraphics(fig3, fullfile(outDir, 'study3_rtd.png'));
