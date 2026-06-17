from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from math import atan2, cos, exp, pi, sin, sqrt
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from induction_motor_model import (
    MotorParams,
    StationaryAlphaBetaInductionMotor,
    abc_to_alphabeta,
    alphabeta_to_abc,
)
from load_model import rpm_to_rad_per_sec
from vector_controller import wrap_angle


ROOT = Path(__file__).resolve().parent
DEFAULT_PLOT = ROOT / "continuous_foc_result.png"
DEFAULT_CSV = ROOT / "continuous_foc_result.csv"


@dataclass
class ContinuousFocConfig:
    torque_ref: float = -220.0
    speed_rpm: float = 5000.0
    stop_time: float = 0.12
    dt: float = 1e-6
    control_period: float = 1e-6
    phase_current_max: float = 550.0
    torque_max: float = 220.0
    current_omega_cc: float = 1000.0
    current_kp_scale: float = 1.0
    current_ki_scale: float = 1.0
    torque_filter_tau: float = 0.01
    flux_time_constant: float = 0.0
    flux_calc_input: str = "command"
    slip_iq: str = "command"
    slip_flux: str = "feedback"
    decoupling_current: str = "command"
    decoupling_flux: str = "feedback"
    decoupling_gain: float = 1.0
    voltage_limit: float = 0.0
    use_abs_flux_for_slip: bool = True
    use_anti_windup: bool = True


class PI:
    def __init__(self) -> None:
        self.integral = 0.0

    def step(self, error: float, kp: float, ki: float, dt: float, limit: float, anti_windup: bool) -> float:
        cand_i = self.integral + error * dt
        cand = kp * (error + ki * cand_i)
        if limit > 0.0:
            out = max(-limit, min(limit, cand))
            if anti_windup:
                saturated = cand != out
                drives_deeper = (cand > limit and error > 0.0) or (cand < -limit and error < 0.0)
                if not (saturated and drives_deeper):
                    self.integral = cand_i
            else:
                self.integral = cand_i
            return out
        self.integral = cand_i
        return kp * (error + ki * self.integral)


class DirectFoc:
    def __init__(self, motor_params: MotorParams, cfg: ContinuousFocConfig) -> None:
        self.p = motor_params
        self.cfg = cfg
        self.id_pi = PI()
        self.iq_pi = PI()
        self.theta_e = 0.0
        self.slip_angle = 0.0
        self.rotor_flux_hat = 0.0
        self.filtered_torque_ref = 0.0
        self.last = {}

    def current_gains(self) -> tuple[float, float]:
        sigma = 1.0 - (self.p.lm * self.p.lm) / max(1e-12, self.p.ls * self.p.lr)
        sigma_ls = max(1e-12, sigma * self.p.ls)
        kp = sigma_ls * self.cfg.current_omega_cc * self.cfg.current_kp_scale
        ki = self.p.rs / sigma_ls * self.cfg.current_ki_scale
        return kp, ki

    def refs(self) -> tuple[float, float, float, float]:
        if self.cfg.torque_filter_tau > 0.0:
            a = 1.0 - exp(-self.cfg.control_period / self.cfg.torque_filter_tau)
            self.filtered_torque_ref += a * (self.cfg.torque_ref - self.filtered_torque_ref)
            torque_ref = self.filtered_torque_ref
        else:
            torque_ref = self.cfg.torque_ref
            self.filtered_torque_ref = torque_ref
        ratio = max(0.0, min(1.0, abs(torque_ref) / max(1e-9, abs(self.cfg.torque_max))))
        axis_current = max(0.0, self.cfg.phase_current_max) * sqrt(3.0) * ratio / sqrt(2.0)
        id_ref = axis_current
        iq_ref = axis_current if torque_ref >= 0.0 else -axis_current
        return id_ref, iq_ref, torque_ref, self.p.lm * id_ref

    def step(self, ia: float, ib: float, ic: float, omega_m: float) -> tuple[float, float]:
        i_alpha, i_beta = abc_to_alphabeta(ia, ib, ic)
        sin_t = sin(self.theta_e)
        cos_t = cos(self.theta_e)
        id_meas = cos_t * i_alpha + sin_t * i_beta
        iq_meas = -sin_t * i_alpha + cos_t * i_beta

        id_ref, iq_ref, torque_ref, flux_ref = self.refs()
        tau_flux = self.cfg.flux_time_constant if self.cfg.flux_time_constant > 0.0 else self.p.lr / max(1e-9, self.p.rr)
        flux_id = id_ref if self.cfg.flux_calc_input == "command" else id_meas
        self.rotor_flux_hat += ((self.p.lm * flux_id - self.rotor_flux_hat) / tau_flux) * self.cfg.control_period

        flux_for_slip = flux_ref if self.cfg.slip_flux == "command" else self.rotor_flux_hat
        if self.cfg.use_abs_flux_for_slip:
            slip_den_flux = abs(flux_for_slip)
        else:
            slip_den_flux = flux_for_slip
        if abs(slip_den_flux) < 1e-9:
            slip_den_flux = 1e-9 if slip_den_flux >= 0.0 else -1e-9
        iq_for_slip = iq_ref if self.cfg.slip_iq == "command" else iq_meas
        slip_omega = self.p.rr * self.p.lm * iq_for_slip / (self.p.lr * slip_den_flux)
        omega_e = self.p.pole_pairs * omega_m + slip_omega

        kp, ki = self.current_gains()
        limit = self.cfg.voltage_limit
        vd_pi = self.id_pi.step(id_ref - id_meas, kp, ki, self.cfg.control_period, limit, self.cfg.use_anti_windup)
        vq_pi = self.iq_pi.step(iq_ref - iq_meas, kp, ki, self.cfg.control_period, limit, self.cfg.use_anti_windup)

        sigma = 1.0 - (self.p.lm * self.p.lm) / max(1e-12, self.p.ls * self.p.lr)
        sigma_ls = sigma * self.p.ls
        id_ff = id_ref if self.cfg.decoupling_current == "command" else id_meas
        iq_ff = iq_ref if self.cfg.decoupling_current == "command" else iq_meas
        flux_ff = flux_ref if self.cfg.decoupling_flux == "command" else self.rotor_flux_hat
        vd_ff = -omega_e * sigma_ls * iq_ff
        vq_ff = omega_e * (sigma_ls * id_ff + (self.p.lm / self.p.lr) * flux_ff)
        vd = vd_pi + self.cfg.decoupling_gain * vd_ff
        vq = vq_pi + self.cfg.decoupling_gain * vq_ff
        if limit > 0.0:
            mag = sqrt(vd * vd + vq * vq)
            if mag > limit:
                scale = limit / mag
                vd *= scale
                vq *= scale

        self.slip_angle = wrap_angle(self.slip_angle + slip_omega * self.cfg.control_period)
        self.theta_e = wrap_angle(self.theta_e + omega_e * self.cfg.control_period)

        v_alpha = cos_t * vd - sin_t * vq
        v_beta = sin_t * vd + cos_t * vq
        self.last = {
            "id": id_meas,
            "iq": iq_meas,
            "id_ref": id_ref,
            "iq_ref": iq_ref,
            "torque_ref": torque_ref,
            "flux_ref": flux_ref,
            "flux_hat": self.rotor_flux_hat,
            "slip_omega": slip_omega,
            "omega_e": omega_e,
            "theta_e": self.theta_e,
            "vd": vd,
            "vq": vq,
            "vd_pi": vd_pi,
            "vq_pi": vq_pi,
            "vd_ff": self.cfg.decoupling_gain * vd_ff,
            "vq_ff": self.cfg.decoupling_gain * vq_ff,
        }
        return v_alpha, v_beta


def flux_angle(motor: StationaryAlphaBetaInductionMotor) -> float:
    return atan2(motor.state.psi_rq, motor.state.psi_rd)


def dq(alpha: float, beta: float, theta: float) -> tuple[float, float]:
    return cos(theta) * alpha + sin(theta) * beta, -sin(theta) * alpha + cos(theta) * beta


def run_case(cfg: ContinuousFocConfig, name: str) -> list[dict]:
    p = MotorParams()
    motor = StationaryAlphaBetaInductionMotor(p)
    foc = DirectFoc(p, cfg)
    omega_m = rpm_to_rad_per_sec(cfg.speed_rpm)
    v_alpha = 0.0
    v_beta = 0.0
    next_control = 0.0
    rows: list[dict] = []
    steps = int(cfg.stop_time / cfg.dt)
    sample_stride = max(1, int(50e-6 / cfg.dt))

    for k in range(steps):
        t = k * cfg.dt
        out_now = motor.output(foc.theta_e)
        if t + 0.5 * cfg.dt >= next_control:
            v_alpha, v_beta = foc.step(out_now.ia, out_now.ib, out_now.ic, omega_m)
            while next_control <= t + 0.5 * cfg.dt:
                next_control += cfg.control_period

        out = motor.step(v_alpha, v_beta, omega_m, cfg.dt, foc.theta_e)
        if k % sample_stride == 0:
            theta_flux = flux_angle(motor)
            i_alpha, i_beta = abc_to_alphabeta(out.ia, out.ib, out.ic)
            id_actual, iq_actual = dq(i_alpha, i_beta, theta_flux)
            vd_actual, vq_actual = dq(v_alpha, v_beta, theta_flux)
            delta = wrap_angle(foc.last.get("theta_e", 0.0) - theta_flux)
            phi_r_abs = sqrt(motor.state.psi_rd * motor.state.psi_rd + motor.state.psi_rq * motor.state.psi_rq)
            phi_dr_ctrl, phi_qr_ctrl = dq(motor.state.psi_rd, motor.state.psi_rq, foc.last.get("theta_e", 0.0))
            row = {
                "case": name,
                "t": t,
                "theta_err_deg": delta * 180.0 / pi,
                "torque": out.torque,
                "torque_ref": foc.last.get("torque_ref", 0.0),
                "id": foc.last.get("id", 0.0),
                "iq": foc.last.get("iq", 0.0),
                "id_actual": id_actual,
                "iq_actual": iq_actual,
                "id_ref": foc.last.get("id_ref", 0.0),
                "iq_ref": foc.last.get("iq_ref", 0.0),
                "vd": foc.last.get("vd", 0.0),
                "vq": foc.last.get("vq", 0.0),
                "vd_actual": vd_actual,
                "vq_actual": vq_actual,
                "vd_pi": foc.last.get("vd_pi", 0.0),
                "vq_pi": foc.last.get("vq_pi", 0.0),
                "vd_ff": foc.last.get("vd_ff", 0.0),
                "vq_ff": foc.last.get("vq_ff", 0.0),
                "slip_omega": foc.last.get("slip_omega", 0.0),
                "flux_hat": foc.last.get("flux_hat", 0.0),
                "phi_r_abs": phi_r_abs,
                "phi_dr_ctrl": phi_dr_ctrl,
                "phi_qr_ctrl": phi_qr_ctrl,
            }
            rows.append(row)
    return rows


def plot(rows_by_case: dict[str, list[dict]], path: Path) -> None:
    fig, axes = plt.subplots(7, 2, figsize=(16, 13.5), sharex="col")
    fig.suptitle("Continuous FOC simulator: no inverter, no delay, direct voltage", fontsize=15, fontweight="bold")
    for col, (name, rows) in enumerate(rows_by_case.items()):
        t = [r["t"] * 1000.0 for r in rows]
        axes[0, col].plot(t, [r["theta_err_deg"] for r in rows])
        axes[0, col].set_title(name)
        axes[0, col].set_ylabel("axis error [deg]")
        axes[1, col].plot(t, [r["id_ref"] for r in rows], label="id ref", linestyle="--")
        axes[1, col].plot(t, [r["id"] for r in rows], label="id ctrl")
        axes[1, col].plot(t, [r["id_actual"] for r in rows], label="id flux-axis")
        axes[1, col].set_ylabel("d current [A]")
        axes[2, col].plot(t, [r["iq_ref"] for r in rows], label="iq ref", linestyle="--")
        axes[2, col].plot(t, [r["iq"] for r in rows], label="iq ctrl")
        axes[2, col].plot(t, [r["iq_actual"] for r in rows], label="iq flux-axis")
        axes[2, col].set_ylabel("q current [A]")
        axes[3, col].plot(t, [r["vd"] for r in rows], label="vd")
        axes[3, col].plot(t, [r["vd_actual"] for r in rows], label="vd flux-axis")
        axes[3, col].plot(t, [r["vd_pi"] for r in rows], label="vd PI", linestyle=":")
        axes[3, col].set_ylabel("d voltage [V]")
        axes[4, col].plot(t, [r["vq"] for r in rows], label="vq")
        axes[4, col].plot(t, [r["vq_actual"] for r in rows], label="vq flux-axis")
        axes[4, col].plot(t, [r["vq_pi"] for r in rows], label="vq PI", linestyle=":")
        axes[4, col].set_ylabel("q voltage [V]")
        axes[5, col].plot(t, [r["flux_hat"] for r in rows], label="flux hat")
        axes[5, col].plot(t, [r["phi_r_abs"] for r in rows], label="|phi r|")
        axes[5, col].plot(t, [r["phi_dr_ctrl"] for r in rows], label="phi dr ctrl")
        axes[5, col].plot(t, [r["phi_qr_ctrl"] for r in rows], label="phi qr ctrl")
        axes[5, col].set_ylabel("rotor flux [Wb]")
        axes[6, col].plot(t, [r["torque_ref"] for r in rows], label="torque ref", linestyle="--")
        axes[6, col].plot(t, [r["torque"] for r in rows], label="torque")
        axes[6, col].set_ylabel("torque [Nm]")
        axes[6, col].set_xlabel("time [ms]")
        for ax in axes[:, col]:
            ax.xaxis.set_major_locator(MultipleLocator(25))
            ax.xaxis.set_minor_locator(MultipleLocator(5))
            ax.grid(True, which="major", color="#cbd5e1", alpha=0.8)
            ax.grid(True, which="minor", color="#e2e8f0", alpha=0.7)
            ax.minorticks_on()
            ax.legend(fontsize=7, ncol=3)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=160)


def write_csv(rows_by_case: dict[str, list[dict]], path: Path) -> None:
    rows = [r for case_rows in rows_by_case.values() for r in case_rows]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Continuous induction motor FOC simulator without inverter or delay.")
    parser.add_argument("--stop-time", type=float, default=0.12)
    parser.add_argument("--control-period-us", type=float, default=1.0)
    parser.add_argument("--speed-rpm", type=float, default=5000.0)
    parser.add_argument("--torque", type=float, default=None, help="Run one case only. If omitted, runs +220 and -220 Nm.")
    parser.add_argument("--no-anti-windup", action="store_true")
    parser.add_argument("--no-abs-flux-slip", action="store_true")
    parser.add_argument("--voltage-limit", type=float, default=0.0, help="0 disables voltage vector limit.")
    parser.add_argument("--flux-tau", type=float, default=0.0)
    parser.add_argument("--out", type=Path, default=DEFAULT_PLOT)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    torques = [args.torque] if args.torque is not None else [220.0, -220.0]
    rows_by_case = {}
    for torque in torques:
        cfg = ContinuousFocConfig(
            torque_ref=torque,
            speed_rpm=args.speed_rpm,
            stop_time=args.stop_time,
            control_period=max(1e-6, args.control_period_us * 1e-6),
            voltage_limit=args.voltage_limit,
            flux_time_constant=args.flux_tau,
            use_abs_flux_for_slip=not args.no_abs_flux_slip,
            use_anti_windup=not args.no_anti_windup,
        )
        name = "motoring +220Nm" if torque >= 0 else "regen -220Nm"
        rows_by_case[name] = run_case(cfg, name)
    plot(rows_by_case, args.out)
    write_csv(rows_by_case, args.csv)
    print(args.out)
    print(args.csv)


if __name__ == "__main__":
    main()
