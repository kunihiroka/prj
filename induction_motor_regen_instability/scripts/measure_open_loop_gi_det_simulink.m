%% measure_open_loop_gi_det_simulink.m
% Open-loop frequency-response measurement for an unstable FOC operating point.
%
% Purpose:
%   Measure Gi(jw): [d_vd_PI, d_vq_PI] -> [id_meas, iq_meas]
%   with the current PI loop opened, then compute
%
%       Li(jw) = Ci(jw) * Gi(jw)
%       D(jw)  = det(I + Li(jw))
%
% Required Simulink model setup:
%   1. Hold the operating point externally:
%        speed = 5000 r/min fixed
%        torque/current reference = -220 Nm equivalent
%        theta_ctrl generated from the nominal slip command
%
%   2. Open the current PI loop:
%        disable/ignore current PI outputs
%        feed vd_PI = vd_PI0 + inj(:,1)
%        feed vq_PI = vq_PI0 + inj(:,2)
%
%   3. Add a From Workspace block named/connected to variable:
%        inj_ts
%      where inj_ts has two columns:
%        column 1: d-axis PI-voltage injection [V]
%        column 2: q-axis PI-voltage injection [V]
%
%   4. Log measured dq current to logsout as:
%        id_meas
%        iq_meas
%
% Notes:
%   - Use small injection amplitude to stay in small-signal range.
%   - For unstable closed-loop points, keep the current loop open while measuring Gi.
%   - Injection is delayed by injStartTime so the model can settle before
%     the sine perturbation starts.
%   - This script assumes continuous or fixed-step simulation is already configured.

clear; clc;

%% ===================== User settings =====================

model = "YOUR_SIMULINK_MODEL_NAME";  % <-- change this

% Frequency points. Put dense points around the suspected 3 Hz mode.
f_Hz = unique([ ...
    logspace(log10(0.5), log10(20), 45), ...
    linspace(2.0, 5.0, 31), ...
    logspace(log10(20), log10(1000), 45) ...
]);
f_Hz = sort(f_Hz(:).');
w = 2*pi*f_Hz;

injAmp_V = 0.5;          % small voltage perturbation amplitude [V]
injStartTime = 0.3;      % sine injection start time [s]
settleCycles = 4;        % cycles discarded before fitting
fitCycles = 6;           % cycles used for sine fitting
pointsPerCycle = 200;    % only for injection waveform generation

% Current PI gains used after Gi measurement.
Kp = 0.0751389;          % replace with actual Kp
Ki = 101.412;            % replace with actual Ki [1/s] for Kp*(1+Ki/s)

% Model variables that switch to open-loop measurement mode.
% These variables must be consumed by your Simulink model.
openLoopEnable = true;
vdPI0 = 0.0;             % steady-state PI d-voltage around operating point [V]
vqPI0 = 0.0;             % steady-state PI q-voltage around operating point [V]

assignin("base", "openLoopEnable", openLoopEnable);
assignin("base", "vdPI0", vdPI0);
assignin("base", "vqPI0", vqPI0);

%% ===================== Run measurements =====================

load_system(model);

Gi = zeros(2, 2, numel(f_Hz));
coh = zeros(2, 2, numel(f_Hz));  % simple fit quality indicator

for k = 1:numel(f_Hz)
    f = f_Hz(k);
    omega = 2*pi*f;
    T = 1/f;
    tEnd = injStartTime + (settleCycles + fitCycles) * T;
    dt = T / pointsPerCycle;
    t = (0:dt:tEnd).';

    fprintf("Measuring %.4g Hz (%d/%d)\n", f, k, numel(f_Hz));

    for inputAxis = 1:2
        inj = zeros(numel(t), 2);
        injActive = t >= injStartTime;
        inj(injActive, inputAxis) = injAmp_V * sin(omega*(t(injActive) - injStartTime));
        inj_ts = timeseries(inj, t);
        assignin("base", "inj_ts", inj_ts);

        simIn = Simulink.SimulationInput(model);
        simIn = simIn.setVariable("inj_ts", inj_ts);
        simIn = simIn.setVariable("openLoopEnable", openLoopEnable);
        simIn = simIn.setVariable("vdPI0", vdPI0);
        simIn = simIn.setVariable("vqPI0", vqPI0);
        simIn = simIn.setModelParameter("StopTime", num2str(tEnd));

        out = sim(simIn);

        [tout, id, iq] = extract_logged_currents(out);

        % Fit only final cycles to remove startup transient.
        tFit0 = injStartTime + settleCycles * T;
        idx = tout >= tFit0;
        tt = tout(idx);
        yy = [id(idx), iq(idx)];

        uFit = injAmp_V * sin(omega*(tt - injStartTime));
        U = fit_complex_sine(tt, uFit, omega);
        for outputAxis = 1:2
            Y = fit_complex_sine(tt, yy(:, outputAxis), omega);
            Gi(outputAxis, inputAxis, k) = Y / U;
            coh(outputAxis, inputAxis, k) = fit_quality(tt, yy(:, outputAxis), omega);
        end
    end
end

%% ===================== Compute Li and det(I+Li) =====================

Dchar = zeros(1, numel(f_Hz));
Li = zeros(2, 2, numel(f_Hz));
I2 = eye(2);

for k = 1:numel(f_Hz)
    Ci = Kp * (1 + Ki/(1j*w(k)));
    Li(:,:,k) = Ci * Gi(:,:,k);
    Dchar(k) = det(I2 + Li(:,:,k));
end

%% ===================== Plot =====================

figure("Name", "Open-loop measured Gi and det(I+Li)", "Color", "w");

subplot(3,2,1);
semilogx(f_Hz, squeeze(20*log10(abs(Gi(1,1,:)))), "LineWidth", 1.3); hold on;
semilogx(f_Hz, squeeze(20*log10(abs(Gi(1,2,:)))), "LineWidth", 1.3);
semilogx(f_Hz, squeeze(20*log10(abs(Gi(2,1,:)))), "LineWidth", 1.3);
semilogx(f_Hz, squeeze(20*log10(abs(Gi(2,2,:)))), "LineWidth", 1.3);
grid on; ylabel("|Gi| [dB]"); title("Measured Gi gain");
legend("id/vd", "id/vq", "iq/vd", "iq/vq", "Location", "best");

subplot(3,2,3);
semilogx(f_Hz, squeeze(rad2deg(unwrap(angle(Gi(1,1,:))))), "LineWidth", 1.3); hold on;
semilogx(f_Hz, squeeze(rad2deg(unwrap(angle(Gi(1,2,:))))), "LineWidth", 1.3);
semilogx(f_Hz, squeeze(rad2deg(unwrap(angle(Gi(2,1,:))))), "LineWidth", 1.3);
semilogx(f_Hz, squeeze(rad2deg(unwrap(angle(Gi(2,2,:))))), "LineWidth", 1.3);
grid on; ylabel("phase [deg]"); xlabel("frequency [Hz]"); title("Measured Gi phase");

subplot(3,2,2);
semilogx(f_Hz, 20*log10(abs(Dchar)), "LineWidth", 1.6);
grid on; ylabel("|det(I+Li)| [dB]"); title("Characteristic function gain");
yline(0, "k--");

subplot(3,2,4);
semilogx(f_Hz, rad2deg(unwrap(angle(Dchar))), "LineWidth", 1.6);
grid on; ylabel("phase [deg]"); xlabel("frequency [Hz]"); title("Characteristic function phase");

subplot(3,2,6);
plot(real(Dchar), imag(Dchar), "LineWidth", 1.6); hold on;
plot(0, 0, "rx", "MarkerSize", 10, "LineWidth", 2);
grid on; axis equal;
xlabel("Re{det(I+Li)}"); ylabel("Im{det(I+Li)}");
title("Nyquist of det(I+Li); danger point = origin");

subplot(3,2,5);
semilogx(f_Hz, squeeze(min(coh, [], [1 2])), "LineWidth", 1.3);
grid on; ylim([0 1.05]);
ylabel("min fit quality"); xlabel("frequency [Hz]");
title("Sine-fit quality check");

%% ===================== Local functions =====================

function [t, id, iq] = extract_logged_currents(simOut)
    logsout = simOut.logsout;
    idSig = logsout.get("id_meas");
    iqSig = logsout.get("iq_meas");
    t = idSig.Values.Time(:);
    id = idSig.Values.Data(:);
    iq = iqSig.Values.Data(:);
end

function X = fit_complex_sine(t, x, omega)
    % Fit x(t) = a*sin(wt) + b*cos(wt) + c.
    % Complex phasor is X = b - j*a for convention Re{X exp(jwt)}.
    t = t(:);
    x = x(:);
    H = [sin(omega*t), cos(omega*t), ones(size(t))];
    theta = H \ x;
    a = theta(1);
    b = theta(2);
    X = b - 1j*a;
end

function q = fit_quality(t, x, omega)
    % R^2-like quality for single-frequency fit.
    t = t(:);
    x = x(:);
    H = [sin(omega*t), cos(omega*t), ones(size(t))];
    theta = H \ x;
    xhat = H * theta;
    ssRes = sum((x - xhat).^2);
    ssTot = sum((x - mean(x)).^2) + eps;
    q = max(0, min(1, 1 - ssRes/ssTot));
end
