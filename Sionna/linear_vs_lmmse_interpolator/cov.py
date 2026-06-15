"""Visualize TDL-A frequency/time covariance matrices vs delay spread and UE speed."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
from sionna.phy.nr import PUSCHConfig, PUSCHTransmitter
from sionna.phy.ofdm import tdl_freq_cov_mat, tdl_time_cov_mat

NUM_LEVELS = 9

# 1 PRB in frequency x 1 slot in time (matches typical NR DMRS slot)
NUM_SUBCARRIERS = 12
NUM_OFDM_SYMBOLS = 14

# NR timing from PUSCH grid (15 kHz SCS, 14-symbol slot)
rg = PUSCHTransmitter(PUSCHConfig()).resource_grid
subcarrier_spacing = rg.subcarrier_spacing
ofdm_symbol_duration = rg.ofdm_symbol_duration

carrier_frequency = 3.5e9
SPEED_OF_LIGHT = 299_792_458.0
delay_spreads_ns = np.round(np.linspace(0, 2000, NUM_LEVELS)).astype(int).tolist()
speeds_ms = np.round(np.linspace(0, 50, NUM_LEVELS), 1).tolist()

DOPPLER_NOTE = (
    "Note: max Doppler f_D = v*f_c/c is shown in speed labels; "
    f"f_c = {carrier_frequency / 1e9:.1f} GHz."
)


def doppler_hz(speed_ms: float) -> float:
    """Max one-sided Doppler [Hz]: f_D = v * f_c / c."""
    return speed_ms * carrier_frequency / SPEED_OF_LIGHT


def speed_with_doppler_label(speed_ms: float) -> str:
    f_d = doppler_hz(speed_ms)
    return f"{speed_ms:g} m/s ({f_d:.0f} Hz Doppler)"


def to_numpy(cov: torch.Tensor) -> np.ndarray:
    return cov.detach().cpu().numpy()


def cov_summary_line(name: str, cov: torch.Tensor, param_label: str, param_value: float) -> str:
    c = to_numpy(cov)
    mag = np.abs(c)
    return (
        f"{name:12s} | {param_label}={param_value:10.4g} | "
        f"shape={tuple(cov.shape)} | "
        f"|R|_max={mag.max():.4f} | diag mean={np.mean(np.real(np.diag(c))):.4f}"
    )


def print_cov_summary(name: str, cov: torch.Tensor, param_label: str, param_value: float):
    print(cov_summary_line(name, cov, param_label, param_value))


def _fmt_complex(z: complex, width: int = 20) -> str:
    if abs(z.imag) < 1e-12:
        return f"{z.real: .8f}".rjust(width)
    sign = "+" if z.imag >= 0 else "-"
    return f"{z.real: .6f}{sign}{abs(z.imag):.6f}j".rjust(width)


def format_cov_matrix(cov: torch.Tensor, label: str) -> str:
    """Format covariance as an n x m complex matrix grid."""
    c = to_numpy(cov)
    n, m = c.shape
    col_w = 20
    lines = [f"{label} ({n} x {m} complex):"]

    header = "      " + "".join(f"[{j:2d}]".center(col_w) for j in range(m))
    lines.append(header)
    lines.append("      " + "-" * (col_w * m))

    for i in range(n):
        row = f"[{i:2d}]  "
        row += "".join(_fmt_complex(c[i, j], col_w) for j in range(m))
        lines.append(row)

    return "\n".join(lines)


def write_cov_txt(path: Path, freq_mats, time_mats, freq_params, time_params, speed_labels):
    with path.open("w", encoding="utf-8") as f:
        f.write("TDL-A covariance matrices\n")
        f.write(f"num_levels = {NUM_LEVELS}\n")
        f.write(f"num_subcarriers = {NUM_SUBCARRIERS}\n")
        f.write(f"num_ofdm_symbols = {NUM_OFDM_SYMBOLS}\n")
        f.write(f"subcarrier_spacing_Hz = {subcarrier_spacing}\n")
        f.write(f"ofdm_symbol_duration_s = {ofdm_symbol_duration}\n")
        f.write(f"carrier_frequency_Hz = {carrier_frequency}\n")
        f.write(f"{DOPPLER_NOTE}\n\n")

        f.write(f"=== Frequency covariance R^(f) — {NUM_SUBCARRIERS} x {NUM_SUBCARRIERS} complex (delay spread) ===\n")
        for ds_ns, cov in zip(freq_params, freq_mats):
            f.write(f"\n--- delay_spread_ns = {ds_ns} ---\n")
            f.write(cov_summary_line("R^(f)", cov, "delay_spread_ns", ds_ns) + "\n\n")
            f.write(format_cov_matrix(cov, "R^(f)") + "\n")

        f.write(
            f"\n=== Time covariance R^(t) — {NUM_OFDM_SYMBOLS} x {NUM_OFDM_SYMBOLS} "
            f"complex (UE speed) ===\n"
        )
        for speed, label, cov in zip(time_params, speed_labels, time_mats):
            f.write(f"\n--- {label} ---\n")
            f.write(cov_summary_line("R^(t)", cov, "speed_m/s", speed) + "\n")
            f.write(f"             | f_D = {doppler_hz(speed):.2f} Hz\n\n")
            f.write(format_cov_matrix(cov, "R^(t)") + "\n")


def plot_cov_grid(matrices, param_values, param_label, title, xlabel, ylabel, note=None):
    n = len(matrices)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.2 * ncols, 3.8 * nrows), squeeze=False, layout="constrained"
    )

    vmax = max(np.abs(to_numpy(m)).max() for m in matrices)
    im = None

    for idx, (cov, pval) in enumerate(zip(matrices, param_values)):
        ax = axes[idx // ncols, idx % ncols]
        c = to_numpy(cov)
        im = ax.imshow(np.abs(c), origin="upper", aspect="auto", vmin=0, vmax=vmax, cmap="viridis")
        ax.set_title(f"{param_label} = {pval}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    for idx in range(n, nrows * ncols):
        axes[idx // ncols, idx % ncols].axis("off")

    if im is not None:
        fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85, label="|R|")
    fig.suptitle(title)
    if note:
        fig.supxlabel(note, fontsize=9)
    return fig


# --- Frequency covariance vs delay spread (fixed geometry) ---
print(f"=== Frequency covariance R^(f) [{NUM_SUBCARRIERS} x {NUM_SUBCARRIERS}] vs delay spread ===")
print(
    f"subcarrier_spacing = {subcarrier_spacing/1e3:.0f} kHz, "
    f"fft_size = {NUM_SUBCARRIERS} (1 PRB), levels = {NUM_LEVELS}"
)
freq_mats = []
for ds_ns in delay_spreads_ns:
    ds_s = ds_ns * 1e-9
    cov = tdl_freq_cov_mat(
        "A", subcarrier_spacing, NUM_SUBCARRIERS, ds_s
    ).to(torch.complex64)
    freq_mats.append(cov)
    print_cov_summary("R^(f)", cov, "delay_spread_ns", ds_ns)

# --- Time covariance vs UE speed (fixed geometry) ---
print(f"\n=== Time covariance R^(t) [{NUM_OFDM_SYMBOLS} x {NUM_OFDM_SYMBOLS}] vs UE speed ===")
print(
    f"carrier_frequency = {carrier_frequency/1e9:.1f} GHz, "
    f"T_sym = {ofdm_symbol_duration*1e6:.2f} us, "
    f"num_ofdm_symbols = {NUM_OFDM_SYMBOLS}, levels = {NUM_LEVELS}"
)
print(DOPPLER_NOTE)
time_mats = []
speed_labels = [speed_with_doppler_label(s) for s in speeds_ms]
for speed, label in zip(speeds_ms, speed_labels):
    cov = tdl_time_cov_mat(
        "A", speed, carrier_frequency, ofdm_symbol_duration, NUM_OFDM_SYMBOLS
    ).to(torch.complex64)
    time_mats.append(cov)
    print_cov_summary("R^(t)", cov, "speed_m/s", speed)
    print(f"             | {label}")

out_dir = Path(__file__).resolve().parent
cov_txt_path = out_dir / "cov.txt"
write_cov_txt(cov_txt_path, freq_mats, time_mats, delay_spreads_ns, speeds_ms, speed_labels)
print(f"\nSaved covariance data to {cov_txt_path}")

img_dir = out_dir / "imgs"
img_dir.mkdir(exist_ok=True)

fig_freq = plot_cov_grid(
    freq_mats,
    [f"{v} ns" for v in delay_spreads_ns],
    "Delay spread",
    f"Frequency covariance |R^(f)| vs delay spread (TDL-A, 12 SC, {NUM_LEVELS} levels)",
    "Subcarrier index",
    "Subcarrier index",
)
fig_freq.savefig(img_dir / "cov_freq_delay_spread.png", dpi=150)

fig_time = plot_cov_grid(
    time_mats,
    speed_labels,
    "UE speed",
    f"Time covariance |R^(t)| vs UE speed (TDL-A, 14 symbols, {NUM_LEVELS} levels)",
    "OFDM symbol index",
    "OFDM symbol index",
    note=DOPPLER_NOTE,
)
fig_time.savefig(img_dir / "cov_time_speed.png", dpi=150)

print(f"\nSaved figures to {img_dir}/")
plt.show()
