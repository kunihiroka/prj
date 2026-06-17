%% plot_l_delta_bode_matlab.m
% Axis-error -> axis-error open-loop Bode plot for induction-motor FOC.
%
% This script computes the reduced axis-error loop
%
%   L_delta(jw)
%     = - S_omega * Gm(jw) * Ci(jw) * M_delta_e / (jw)
%
% where
%
%   Δδ
%    -> M_delta_e       : axis error to current-control error
%    -> Ci              : current PI
%    -> Gm              : motor voltage to [id, iq, phi_d]
%    -> S_omega         : [Δid, Δiq, Δphi] to real slip variation
%    -> -               : Δδdot = - Δomega_slip_real
%    -> 1/(jw)          : angle integration
%
% Because this definition already includes
%
%   Δδdot = - Δomega_slip_real,
%
% the characteristic equation of this reduced self-coupling model is
%
%   1 - L_delta(s) = 0.
%
% Therefore the danger point is L_delta = +1, or equivalently
% D_delta = 1 - L_delta = 0.
%
% It compares motoring and regeneration at 5000 r/min, 220 Nm equivalent.

clear; clc;

%% Motor and control constants
Rs  = 0.00762;
Rr  = 0.008041;
Lls = 0.0000419;
Llr = 0.0000419;
Lm  = 0.0001583;

Ls = Lls + Lm;
Lr = Llr + Lm;
polePairs = 4;

speedRpm = 5000.0;
omegaM = speedRpm * 2*pi/60;
omegaR = polePairs * omegaM;

sigma = 1 - Lm^2/(Ls*Lr);
sigmaLs = sigma * Ls;

omegaCc = 1000.0;
Kp = sigmaLs * omegaCc;
Ki = Rs / sigmaLs;             % for Ci = Kp*(1 + Ki/s)

IphaseMax = 550.0;
Id0 = IphaseMax * sqrt(3) / sqrt(2);
IqMot = +Id0;
IqReg = -Id0;

Kslip = Rr * Lm / Lr;

%% Frequency vector
f = logspace(-1, 4, 1800);     % Hz
w = 2*pi*f;                    % rad/s

%% Compute L_delta
Lmot = calc_l_delta(IqMot, w, Rs, Rr, Ls, Lr, Lm, omegaR, Kp, Ki, Id0, Kslip);
Lreg = calc_l_delta(IqReg, w, Rs, Rr, Ls, Lr, Lm, omegaR, Kp, Ki, Id0, Kslip);

%% Summaries
[fcMot, phMot] = first_zero_db_crossing(f, Lmot);
[fcReg, phReg] = first_zero_db_crossing(f, Lreg);

[distMot, idxMot] = min(abs(1 - Lmot));
[distReg, idxReg] = min(abs(1 - Lreg));

fprintf("motoring: 0 dB crossing = %.4g Hz, phase = %.2f deg\n", fcMot, phMot);
fprintf("regen:    0 dB crossing = %.4g Hz, phase = %.2f deg\n", fcReg, phReg);
fprintf("motoring: nearest to +1 / min |1-L| at %.4g Hz, |1-L| = %.4g, L = %.4g%+.4gj\n", ...
    f(idxMot), distMot, real(Lmot(idxMot)), imag(Lmot(idxMot)));
fprintf("regen:    nearest to +1 / min |1-L| at %.4g Hz, |1-L| = %.4g, L = %.4g%+.4gj\n", ...
    f(idxReg), distReg, real(Lreg(idxReg)), imag(Lreg(idxReg)));

%% Plot
figure("Name", "Axis-error loop L_delta", "Color", "w");

subplot(2,2,1);
semilogx(f, 20*log10(abs(Lmot)), "LineWidth", 1.8); hold on;
semilogx(f, 20*log10(abs(Lreg)), "--", "LineWidth", 1.8);
yline(0, "k-");
grid on;
xlabel("frequency [Hz]");
ylabel("|L_\delta| [dB]");
title("Bode gain");
legend("motoring +220 Nm", "regen -220 Nm", "Location", "best");

subplot(2,2,3);
semilogx(f, unwrap(angle(Lmot))*180/pi, "LineWidth", 1.8); hold on;
semilogx(f, unwrap(angle(Lreg))*180/pi, "--", "LineWidth", 1.8);
yline(-180, ":", "Color", [0.45 0.5 0.6]);
yline(0, "k-");
grid on;
xlabel("frequency [Hz]");
ylabel("phase [deg]");
title("Bode phase");
legend("motoring +220 Nm", "regen -220 Nm", "Location", "best");

subplot(2,2,2);
plot(real(Lmot), imag(Lmot), "LineWidth", 1.8); hold on;
plot(real(Lreg), imag(Lreg), "--", "LineWidth", 1.8);
plot(1, 0, "rx", "MarkerSize", 10, "LineWidth", 2);
grid on; axis equal;
xlim([-1.0 2.5]);
ylim([-1.75 1.75]);
xlabel("Re{L_\delta(j\omega)}");
ylabel("Im{L_\delta(j\omega)}");
title("Nyquist of L_\delta; danger point = +1");
legend("motoring", "regen", "+1+j0", "Location", "best");

subplot(2,2,4);
plot(real(1 - Lmot), imag(1 - Lmot), "LineWidth", 1.8); hold on;
plot(real(1 - Lreg), imag(1 - Lreg), "--", "LineWidth", 1.8);
plot(0, 0, "rx", "MarkerSize", 10, "LineWidth", 2);
grid on; axis equal;
xlim([-1.5 3.0]);
ylim([-1.75 1.75]);
xlabel("Re{1-L_\delta(j\omega)}");
ylabel("Im{1-L_\delta(j\omega)}");
title("Characteristic D_\delta = 1-L_\delta");
legend("motoring", "regen", "0+j0", "Location", "best");

%% ===================== Local functions =====================

function Ldelta = calc_l_delta(iq0, w, Rs, Rr, Ls, Lr, Lm, omegaR, Kp, Ki, Id0, Kslip)
    [A, B, C] = linearize_motor(iq0, Rs, Rr, Ls, Lr, Lm, omegaR, Id0);

    phi = Lm * Id0;
    Somega = [0, Kslip/phi, -Kslip*iq0/(phi^2)];
    MdeltaE = [-iq0; Id0];

    n = size(A, 1);
    I = eye(n);
    Ldelta = zeros(size(w));

    for k = 1:numel(w)
        jw = 1j*w(k);
        Gm = C * ((jw*I - A) \ B);          % [id, iq, phi] / [vd, vq]
        Ci = Kp * (1 + Ki/jw);              % scalar PI
        Hdot = -Somega * Gm * (Ci * MdeltaE);
        Ldelta(k) = Hdot / jw;              % δdot -> δ integration
    end
end

function [A, B, C] = linearize_motor(iq0, Rs, Rr, Ls, Lr, Lm, omegaR, Id0)
    [x0, slip, omegaE] = operating_state(iq0, Rs, Rr, Ls, Lr, Lm, omegaR, Id0);

    n = 4;
    A = zeros(n, n);
    B = zeros(n, 2);

    for j = 1:n
        h = 1e-7 * max(1, abs(x0(j)));
        xp = x0;
        xm = x0;
        xp(j) = xp(j) + h;
        xm(j) = xm(j) - h;
        A(:, j) = (motor_rhs(xp, [0;0], slip, omegaE, Rs, Rr, Ls, Lr, Lm) ...
                 - motor_rhs(xm, [0;0], slip, omegaE, Rs, Rr, Ls, Lr, Lm)) / (2*h);
    end

    for j = 1:2
        h = 1e-4;
        up = [0;0];
        um = [0;0];
        up(j) = +h;
        um(j) = -h;
        B(:, j) = (motor_rhs(x0, up, slip, omegaE, Rs, Rr, Ls, Lr, Lm) ...
                 - motor_rhs(x0, um, slip, omegaE, Rs, Rr, Ls, Lr, Lm)) / (2*h);
    end

    detL = Ls*Lr - Lm^2;
    C = [
        Lr/detL, 0,       -Lm/detL, 0;       % id
        0,       Lr/detL, 0,        -Lm/detL; % iq
        0,       0,       1,        0         % phi_d
    ];
end

function [x0, slip, omegaE] = operating_state(iq0, Rs, Rr, Ls, Lr, Lm, omegaR, Id0)
    %#ok<INUSD> Rs is unused here but kept for consistent function signature.
    phi = Lm * Id0;

    iRd = 0.0;
    iRq = -Lm/Lr * iq0;

    psiSd = Ls*Id0 + Lm*iRd;
    psiSq = Ls*iq0 + Lm*iRq;
    psiRd = phi;
    psiRq = 0.0;

    slip = Rr*Lm*iq0 / (Lr*phi);
    omegaE = omegaR + slip;

    x0 = [psiSd; psiSq; psiRd; psiRq];
end

function dx = motor_rhs(x, u, slip, omegaE, Rs, Rr, Ls, Lr, Lm)
    [ids, iqs, idr, iqr] = currents_from_flux(x, Ls, Lr, Lm);

    vd = u(1);
    vq = u(2);

    psiSd = x(1);
    psiSq = x(2);
    psiRd = x(3);
    psiRq = x(4);

    dpsiSd = vd - Rs*ids + omegaE*psiSq;
    dpsiSq = vq - Rs*iqs - omegaE*psiSd;
    dpsiRd = -Rr*idr + slip*psiRq;
    dpsiRq = -Rr*iqr - slip*psiRd;

    dx = [dpsiSd; dpsiSq; dpsiRd; dpsiRq];
end

function [ids, iqs, idr, iqr] = currents_from_flux(x, Ls, Lr, Lm)
    psiSd = x(1);
    psiSq = x(2);
    psiRd = x(3);
    psiRq = x(4);

    detL = Ls*Lr - Lm^2;

    ids = ( Lr*psiSd - Lm*psiRd) / detL;
    iqs = ( Lr*psiSq - Lm*psiRq) / detL;
    idr = (-Lm*psiSd + Ls*psiRd) / detL;
    iqr = (-Lm*psiSq + Ls*psiRq) / detL;
end

function [fc, ph] = first_zero_db_crossing(f, L)
    magDb = 20*log10(abs(L));
    phaseDeg = unwrap(angle(L))*180/pi;

    fc = NaN;
    ph = NaN;
    for k = 1:(numel(f)-1)
        if magDb(k) == 0 || magDb(k)*magDb(k+1) < 0
            a = abs(magDb(k));
            b = abs(magDb(k+1));
            r = a/(a+b);
            fc = f(k) + (f(k+1)-f(k))*r;
            ph = phaseDeg(k) + (phaseDeg(k+1)-phaseDeg(k))*r;
            return;
        end
    end
end
