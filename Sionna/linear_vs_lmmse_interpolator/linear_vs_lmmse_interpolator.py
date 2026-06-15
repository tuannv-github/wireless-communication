from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import torch
from sionna.phy.channel import OFDMChannel
from sionna.phy.channel.tr38901 import TDL
from sionna.phy.nr import PUSCHConfig, PUSCHLSChannelEstimator, PUSCHTransmitter
from sionna.phy.ofdm import LMMSEInterpolator, tdl_freq_cov_mat, tdl_time_cov_mat

# NR PUSCH resource grid with DMRS pilots (1 CDM group → 6 DMRS REs/PRB)
pusch_config = PUSCHConfig()
pusch_config.dmrs.additional_position = 1
pusch_config.dmrs.num_cdm_groups_without_data = 1
transmitter = PUSCHTransmitter(pusch_config)
rg = transmitter.resource_grid

carrier_frequency = 3.5e9
delay_spread = 100e-9
batch_size = 16
no = torch.tensor(0.01)

estimator_kwargs = dict(
    resource_grid=rg,
    dmrs_length=pusch_config.dmrs.length,
    dmrs_additional_position=pusch_config.dmrs.additional_position,
    num_cdm_groups_without_data=pusch_config.dmrs.num_cdm_groups_without_data,
)
estimator_lin = PUSCHLSChannelEstimator(**estimator_kwargs, interpolation_type="lin")


def full_grid_mse(h_est, h_true):
    h_est = h_est.squeeze(3).squeeze(3)
    return torch.mean(torch.abs(h_est - h_true) ** 2).item()


def evaluate_mse(delay_spread_s, speed_ms):
    cov_mat_freq = tdl_freq_cov_mat(
        "A", rg.subcarrier_spacing, rg.fft_size, delay_spread_s
    ).to(torch.complex64)
    cov_mat_time = tdl_time_cov_mat(
        "A", speed_ms, carrier_frequency, rg.ofdm_symbol_duration, rg.num_ofdm_symbols
    ).to(torch.complex64)
    estimator_lmmse = PUSCHLSChannelEstimator(
        **estimator_kwargs,
        interpolator=LMMSEInterpolator(
            rg.pilot_pattern, cov_mat_time, cov_mat_freq, order="f-t"
        ),
    )

    tx, _ = transmitter(batch_size)
    tdl = TDL(
        model="A",
        delay_spread=delay_spread_s,
        carrier_frequency=carrier_frequency,
        min_speed=speed_ms,
        max_speed=speed_ms,
    )
    channel = OFDMChannel(channel_model=tdl, resource_grid=rg, return_channel=True)
    y, h_true = channel(tx, no)
    h_true = h_true.squeeze(3).squeeze(3)

    est_lin, _ = estimator_lin(y, no)
    est_lmmse, _ = estimator_lmmse(y, no)
    return full_grid_mse(est_lin, h_true), full_grid_mse(est_lmmse, h_true)


def plot_mse_subplot(ax, x_values, mse_linear, mse_lmmse, xlabel, title):
    ax.plot(x_values, mse_linear, "o-", label="Linear")
    ax.plot(x_values, mse_lmmse, "s-", label="LMMSE")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("MSE (full grid)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()


# Sweep delay spread (fixed speed = 3 m/s)
speed = 3.0
delay_spreads_ns = np.linspace(10, 2000, 20)
mse_linear_ds = []
mse_lmmse_ds = []

for ds in delay_spreads_ns * 1e-9:
    m_lin, m_lmmse = evaluate_mse(ds, speed)
    mse_linear_ds.append(m_lin)
    mse_lmmse_ds.append(m_lmmse)

print("Delay spread [ns] | MSE Linear | MSE LMMSE")
for ds, m_lin, m_lmmse in zip(delay_spreads_ns, mse_linear_ds, mse_lmmse_ds):
    print(f"{ds:14.1f} | {m_lin:10.6f} | {m_lmmse:9.6f}")

# Sweep speed (fixed delay spread = 100 ns)
speeds_ms = np.linspace(0, 30, 20)
mse_linear_spd = []
mse_lmmse_spd = []

for speed_ms in speeds_ms:
    m_lin, m_lmmse = evaluate_mse(delay_spread, speed_ms)
    mse_linear_spd.append(m_lin)
    mse_lmmse_spd.append(m_lmmse)

print("\nSpeed [m/s] | MSE Linear | MSE LMMSE")
for spd, m_lin, m_lmmse in zip(speeds_ms, mse_linear_spd, mse_lmmse_spd):
    print(f"{spd:11.1f} | {m_lin:10.6f} | {m_lmmse:9.6f}")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
plot_mse_subplot(
    axes[0],
    delay_spreads_ns,
    mse_linear_ds,
    mse_lmmse_ds,
    "Delay spread [ns]",
    "MSE vs delay spread (3 m/s)",
)
plot_mse_subplot(
    axes[1],
    speeds_ms,
    mse_linear_spd,
    mse_lmmse_spd,
    "Speed [m/s]",
    "MSE vs speed (100 ns)",
)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.suptitle("Channel estimation MSE: Linear vs LMMSE (TDL-A)", y=0.98)
img_dir = Path(__file__).resolve().parent / "imgs"
img_dir.mkdir(exist_ok=True)
fig.savefig(img_dir / "mse_comparison.png", dpi=150)
plt.show()
