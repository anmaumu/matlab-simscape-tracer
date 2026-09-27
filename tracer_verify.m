%% TRACER_VERIFY  Verification of the extended tracer model (tracer_net)
outDir = fileparts(mfilename('fullpath'));

fluid = struct('rho0', 850, 'pRef', 0.101325e6, 'beta', 1.5e9, ...
    'nu', 3.2e-5, 'Dc', 1e-9);                      % ISO VG32-like oil
orifice = struct('type', 'orifice', 'Cd', 0.64, 'SR', 1e-5, 'S', inf, ...
    'recovery', 0, 'Recr', 150);
pipe = struct('type', 'pipe', 'L', 1, 'D', 0.01, 'rough', 15e-6, 'Leq', 0, ...
    'lam', 64, 'ReL', 2000, 'ReT', 4000);

%% 1. Orifice (Local Restriction IL equations)
fprintf('--- 1. Orifice ---\n');
rho = fluid.rho0;
dp = logspace(0, 7, 400);
orf = @(E, d) tracer_net('orificeFlow', d, rho, E, fluid);
m0 = arrayfun(@(d) orf(orifice, d), dp);

% Re of the turbulent asymptote at dp = dp_cr must equal Re_cr
Dh = sqrt(4 * orifice.SR / pi);
dpCr = pi/4 * rho / (2 * orifice.SR) * (orifice.Recr * fluid.nu / orifice.Cd)^2;
mTurbAtCr = orifice.Cd * orifice.SR * sqrt(2 * rho * dpCr);
ReAtCr = mTurbAtCr / (rho * orifice.SR) * Dh / fluid.nu;
fprintf('dp_cr = %.1f Pa, Re(turbulent asymptote at dp_cr) = %.2f (Re_cr = %g)\n', dpCr, ReAtCr, orifice.Recr);

% Log-log slope: laminar -> 1, turbulent -> 0.5
slope = diff(log(m0)) ./ diff(log(dp));
fprintf('slope d(log m)/d(log dp): at dp = %.0f Pa -> %.3f, at dp = %.0e Pa -> %.3f\n', ...
    dp(1), slope(1), dp(end), slope(end));

% Area ratio and pressure recovery
oRatio = orifice;  oRatio.S = 2 * orifice.SR;      % S_R/S = 0.5
oRec = oRatio;     oRec.recovery = 1;
dpTest = 1e5;
fprintf('at dp = 1e5 Pa, S_R/S = 0.5: m(no ratio) %.5f, m(ratio) %.5f, m(ratio+recovery) %.5f kg/s\n', ...
    orf(orifice, dpTest), orf(oRatio, dpTest), orf(oRec, dpTest));

fig1 = figure('Name', 'Orifice', 'Position', [100 100 700 450]);
loglog(dp, m0, dp, arrayfun(@(d) orf(oRatio, d), dp), dp, arrayfun(@(d) orf(oRec, d), dp), 'LineWidth', 1.2);
xline(dpCr, ':', 'dp_{cr}'); grid on;
xlabel('\Deltap [Pa]'); ylabel('mdot [kg/s]');
legend('S_R/S = 0', 'S_R/S = 0.5', 'S_R/S = 0.5 + recovery', 'Location', 'northwest');
title('Orifice: laminar (slope 1) to turbulent (slope 1/2)');
exportgraphics(fig1, fullfile(outDir, 'verify1_orifice.png'));

%% 2. Pipe friction
fprintf('\n--- 2. Pipe friction ---\n');
S = pi * pipe.D^2 / 4;
Re = logspace(1, 5, 400);
mRe = Re * S * rho * fluid.nu / pipe.D;
dpPipe = arrayfun(@(mm) tracer_net('pipeDp', mm, rho, pipe, fluid), mRe);
fEff = dpPipe * 2 * rho * pipe.D * S^2 ./ (pipe.L * mRe.^2);
fHaa = arrayfun(@(x) 1 / (-1.8 * log10(6.9 / x + (pipe.rough / pipe.D / 3.7)^1.11))^2, Re);
lam = Re < pipe.ReL;
tur = Re > pipe.ReT;
fprintf('laminar   (Re < %d): max |f/(64/Re) - 1| = %.2e\n', pipe.ReL, max(abs(fEff(lam) ./ (64 ./ Re(lam)) - 1)));
fprintf('turbulent (Re > %d): max |f/f_Haaland - 1| = %.2e\n', pipe.ReT, max(abs(fEff(tur) ./ fHaa(tur) - 1)));
fprintf('dp(mdot) strictly increasing: %d\n', all(diff(dpPipe) > 0));

fig2 = figure('Name', 'Pipe friction', 'Position', [100 100 700 450]);
loglog(Re, fEff, 'LineWidth', 1.6); hold on;
loglog(Re(Re < 4000), 64 ./ Re(Re < 4000), '--', Re(Re > 2000), fHaa(Re > 2000), '--', 'LineWidth', 1);
xline([pipe.ReL pipe.ReT], ':'); grid on;
xlabel('Re'); ylabel('Darcy friction factor f');
legend('model', '64/Re', 'Haaland', 'Location', 'northeast');
title('Pipe friction: laminar / smoothed transition / turbulent');
exportgraphics(fig2, fullfile(outDir, 'verify2_pipe.png'));

%% 3. Network: 5 m pipe split into N segments, pulse injection (flow inlet)
fprintf('\n--- 3. Network: pipe L = 5 m, D = 10 mm, flow inlet ---\n');
N = 20;
Lpipe = 5;
seg = pipe;  seg.L = Lpipe / N;
mIn = 0.0356;
P = struct();
P.fluid = fluid;
P.V = S * seg.L * ones(N, 1);
P.elem = [{struct('type', 'source', 'mdot', mIn)}, repmat({seg}, 1, N)];
P.pIn = NaN;  P.cIn = 0;
P.pOut = 0.2e6;  P.cOut = 0;
P.mInj = 0.001;  P.cInj = 1;  P.injVol = 1;
P.tOn = 5;  P.tOff = 6;  P.tEnd = 80;
P.mEps = 1e-6;
r = tracer_net(P);

iPre = find(r.t < P.tOn, 1, 'last');
ReNet = mIn * seg.D / (S * fluid.rho0 * fluid.nu);
dpHP = 32 * fluid.nu * Lpipe * mIn / (seg.D^2 * S);          % Hagen-Poiseuille
dpSim = r.p(iPre, 1) - P.pOut;
fprintf('Re = %.0f, pipe dp: model %.1f Pa, Hagen-Poiseuille %.1f Pa (err %.2e)\n', ...
    ReNet, dpSim, dpHP, dpSim / dpHP - 1);
fprintf('total mass balance error  : %.2e kg (mass in system %.4f kg)\n', r.massBalanceErr, r.massInSystem(1));
fprintf('tracer mass balance error : %.2e kg (injected %.4f kg)\n', r.tracerBalanceErr, r.injectedTracer);

ts = r.t - (P.tOn + P.tOff) / 2;
E = r.mdot_c(:, end) / r.injectedTracer;
tm = trapz(ts, ts .* E);
s2 = trapz(ts, (ts - tm).^2 .* E);
tau = r.massInSystem(iPre) / mIn;
fprintf('RTD: recovery %.4f, mean %.3f s, tau %.3f s, var/tau^2 %.4f (1/N = %.4f)\n', ...
    trapz(ts, E), tm, tau, s2 / tau^2, 1 / N);

fig3 = figure('Name', 'Network', 'Position', [100 100 800 600]);
subplot(2, 1, 1);
plot(r.t, r.c(:, [1 N/2 N]), 'LineWidth', 1.2); grid on;
legend('vol 1', sprintf('vol %d', N/2), sprintf('vol %d', N)); ylabel('c [-]');
title('Tracer pulse travelling through a 5 m pipe (20 segments)');
subplot(2, 1, 2);
plot(r.t, r.mdot(:, [1 2 end]), 'LineWidth', 1.2); grid on;
legend('inlet source', 'pipe seg 1', 'outlet'); ylabel('mdot [kg/s]'); xlabel('t [s]');
exportgraphics(fig3, fullfile(outDir, 'verify3_network.png'));

%% 4. Diffusion at zero flow: two closed volumes linked by a short pipe
fprintf('\n--- 4. Diffusion at zero flow ---\n');
Q = struct();
Q.fluid = fluid;  Q.fluid.Dc = 1e-6;
short = pipe;  short.L = 0.05;
Q.V = [1e-5; 1e-5];
Q.elem = {struct('type', 'source', 'mdot', 0), short, struct('type', 'closed')};
Q.pIn = NaN;  Q.cIn = 0;  Q.pOut = 0.2e6;  Q.cOut = 0;
Q.c0 = [1; 0];
Q.mInj = 0;  Q.cInj = 0;  Q.injVol = 1;
Q.tOn = 0;  Q.tOff = 0;  Q.tEnd = 10000;
Q.mEps = 1e-6;
r4 = tracer_net(Q);

rhoQ = r4.rho(1, 1);
G = rhoQ * Q.fluid.Dc * S / short.L;
tauD = rhoQ * Q.V(1) / (2 * G);
dcTheory = exp(-r4.t / tauD);
dcSim = r4.c(:, 1) - r4.c(:, 2);
fprintf('time constant %.0f s, max |dc_sim - dc_theory| = %.2e, final c = [%.4f %.4f]\n', ...
    tauD, max(abs(dcSim - dcTheory)), r4.c(end, 1), r4.c(end, 2));
fprintf('max |mdot| in pipe = %.2e kg/s, tracer balance error = %.2e kg\n', ...
    max(abs(r4.mdot(:, 2))), r4.tracerBalanceErr);

fig4 = figure('Name', 'Diffusion', 'Position', [100 100 700 400]);
plot(r4.t, r4.c, 'LineWidth', 1.2); hold on;
plot(r4.t, 0.5 + 0.5 * dcTheory, 'k--', r4.t, 0.5 - 0.5 * dcTheory, 'k--');
grid on; xlabel('t [s]'); ylabel('c [-]');
legend('vol 1', 'vol 2', 'theory'); title('Tracer diffusion between two volumes at zero flow');
exportgraphics(fig4, fullfile(outDir, 'verify4_diffusion.png'));
