# LMMSE Channel Interpolation

Notes on the math behind Sionna's [`LMMSEInterpolator`](https://nvlabs.github.io/sionna/phy/api/ofdm/sionna.phy.ofdm.LMMSEInterpolator.html), as used in `linear_vs_lmmse_interpolator.py`.

## Problem setup

After LS estimation at DMRS pilot positions, we have noisy channel samples on a sparse subset of the OFDM resource grid. **Interpolation** fills the grid by predicting channel values at non-pilot REs from pilot REs.

| Method | Idea |
|--------|------|
| **Linear** | Interpolate along time/frequency with fixed weights (ignores channel statistics). |
| **LMMSE** | Wiener weights from **R** and **Σ**; best linear estimator (see [MMSE / LMMSE theory](#mmse-and-lmmse-detailed)) |

Sionna's `LMMSEInterpolator` takes:

- `h_hat` — LS channel estimates at pilots
- `err_var` — LS error variance at each pilot
- `cov_mat_freq` — frequency covariance **R⁽ᶠ⁾** (from `tdl_freq_cov_mat`)
- `cov_mat_time` — time covariance **R⁽ᵗ⁾** (from `tdl_time_cov_mat`)
- `cov_mat_space` — *(optional)* spatial covariance **R⁽ˢ⁾** across receive antennas; required only when `order` includes `"s"`
- `order` — pass order: `"f"`, `"t"`, and optionally `"s"` in any sequence, e.g. `"f-t"` = frequency → scaling → time, `"f-t-s"` = frequency → scaling → time → scaling → spatial

In our script, `order="f-t"` (no spatial pass), `num_cdm_groups_without_data = 1`, and TDL-A covariances depend on **delay spread** and **UE speed**. [Script-specific matrix sizes](#matrix-shapes-from-linear_vs_lmmse_interpolatorpy) and [pass-by-pass formulas](#sionnas-algorithm-orderf-t) follow below.

---

## Where do **R⁽ᶠ⁾** and **R⁽ᵗ⁾** come from?

LMMSE needs **prior knowledge** of channel correlation:

```python
cov_mat_freq = tdl_freq_cov_mat("A", subcarrier_spacing, fft_size, delay_spread)
cov_mat_time = tdl_time_cov_mat("A", speed, carrier_frequency, ofdm_symbol_duration, num_ofdm_symbols)
```

**Why `fft_size`, not “number of subcarriers”?** Sionna’s API names the third argument `fft_size` because **R⁽ᶠ⁾** is indexed by FFT bin \(u,v = 1,\ldots,M\):

\[
R^{(f)}_{u,v} = \sum_{\ell} P_\ell e^{-j 2\pi \tau_\ell \Delta_f (u-v)}
\]

So the matrix is **\(M \times M\)** where \(M\) is the frequency dimension of the channel vector the interpolator works on.

| Concept | In general | In our script (`4` RB PUSCH) |
|---------|------------|------------------------------|
| `fft_size` | OFDM FFT length \(M\) | **48** (`rg.fft_size`) |
| Effective subcarriers | Data/pilot SC actually used | **48** (`rg.num_effective_subcarriers`) |
| Guard / null SC | Often `fft_size >` used SC | **0** (`rg.num_guard_carriers = [0,0]`) |

Here **`fft_size == num_effective_subcarriers == 48`**, so passing `rg.fft_size` is correct and matches the **48-subcarrier row** in pass 1. In a full-cell grid (e.g. 1024-point FFT with guard bands), \(M\) would still be the FFT size that defines subcarrier spacing \(\Delta_f = B/M\); the covariance must match the same indexing as `h_hat` on the resource grid.

The script uses `rg.fft_size` (not a separate “num subcarriers” field) because that is what Sionna’s `ResourceGrid` and `LMMSEInterpolator` use for the frequency axis.

| Covariance | Driven by | Physical meaning |
|------------|-----------|------------------|
| **R⁽ᶠ⁾** | Delay spread $\tau_{\mathrm{rms}}$ | Wider delay spread → **narrower** coherence bandwidth ($B_c \sim 1/\tau_{\mathrm{rms}}$) → **weaker** correlation between subcarriers at a given spacing (more frequency-selective fading). Narrow delay spread → nearly flat channel → all subcarriers almost perfectly correlated. |
| **R⁽ᵗ⁾** | UE speed (Doppler) | Higher speed → channel changes faster across OFDM symbols → weaker time correlation |

If the assumed covariances mismatch the true channel, LMMSE is suboptimal. Linear interpolation ignores them entirely, so it is simpler but often worse when statistics are known.

---

## LMMSE vs linear (intuition for our plots)

**Linear interpolation** treats all gaps equally — like drawing straight lines between pilots.

**LMMSE** weights neighbors using **R⁽ᶠ⁾** / **R⁽ᵗ⁾**:

- **Large delay spread:** coherence bandwidth shrinks → channel varies **more** across frequency → subcarriers are **less** correlated at a given spacing, but **R⁽ᶠ⁾** still captures who is similar to whom; linear interp. uses wrong fixed weights → LMMSE gains are large.
- **High speed:** **R⁽ᵗ⁾** shrinks temporal correlation → harder to predict between DMRS symbols; linear suffers more because it does not adapt.

That is why MSE gaps in `mse_comparison.png` grow with delay spread and mobility.

---

## Important assumptions (from Sionna docs)

1. **Diagonal error covariance** — after each pass, Sionna keeps only $\mathrm{diag}(\cdot)$ of the matrix MSE (uncorrelated errors across subcarriers / symbols / antennas). This is an **approximate** LMMSE when true errors are correlated.
2. **Known covariances** — **R⁽ᶠ⁾**, **R⁽ᵗ⁾** must match the channel model (here TDL-A).
3. **Pass order** — Sionna `order="f-t"` runs frequency → scaling → time. With `"s"`, an extra scaling step and spatial smoothing are added. The string order can be permuted (e.g. `"t-f-s"`).
4. **Numerics** — COD-based solve avoids explicit inversion of ill-conditioned **R**.

---

## Sionna's algorithm (`order="f-t"`)

The resource grid is $\hat{\mathbf{H}} \in \mathbb{C}^{N \times M}$: **N** OFDM symbols, **M** subcarriers. Sionna's `order="f-t"` runs **frequency interpolation → variance scaling → time interpolation**. (With `order="f-t-s"`, two scaling steps and a spatial pass are added.)

### Pass 1 — frequency interpolation (per OFDM symbol)

For each OFDM symbol $n$ that carries pilots ($n \in \{2, 11\}$ in this script):

**Channel estimate**

$$
\hat{\mathbf{h}}_n^{(1)} = \mathbf{A}_n \hat{\mathbf{h}}_n
$$

**Error variance** (diagonal; Sionna keeps only the diagonal of the matrix MSE)

$$
\mathbf{\Sigma}_n^{(1)} = \mathrm{diag}\left( \mathbf{R^{(f)}} - \mathbf{A}_n \mathbf{\Xi}_n \mathbf{R^{(f)}} \right)
$$

Equivalently, for each subcarrier $m$: $[\mathbf{\Sigma}_n^{(1)}]_{m,m} = [\mathbf{R^{(f)}}]_{m,m} - [\mathbf{A}_n \mathbf{\Xi}_n \mathbf{R^{(f)}}]_{m,m}$.

$\mathbf{\Xi}_n$ is an $M \times M$ diagonal matrix that **zeros columns corresponding to subcarriers without pilots** on symbol $n$.

**Interpolation matrix** ($\mathbf{A}_n$ is $M \times M$; solved via COD, not a direct inverse):

$$
\mathbf{A}_n = \bar{\mathbf{A}}_n \mathbf{\Pi}_n^\intercal, \qquad
\bar{\mathbf{A}}_n = \underset{\mathbf{Z} \in \mathbb{C}^{M \times K_n}}{\arg\min} \left\| \mathbf{Z}\left( \mathbf{\Pi}_n^\intercal \mathbf{R^{(f)}} \mathbf{\Pi}_n + \mathbf{\Sigma}_n \right) - \mathbf{R^{(f)}} \mathbf{\Pi}_n \right\|_{\mathrm{F}}^2
$$

| Matrix | Shape | In this script | Meaning |
|--------|-------|----------------|---------|
| $\hat{\mathbf{h}}_n$ (input) | $M \times 1$ | $48 \times 1$ | row on symbol $n$; $K_n=24$ non-zero at pilot SC |
| $\hat{\mathbf{h}}_n^{(1)}$ (output) | $M \times 1$ | $48 \times 1$ | full row after frequency interpolation |
| $\bar{\mathbf{A}}_n$ | $M \times K_n$ | $48 \times 24$ | core interpolation coefficients |
| $\mathbf{A}_n$ | $M \times M$ | $48 \times 48$ | $\bar{\mathbf{A}}_n \mathbf{\Pi}_n^\intercal$ |
| $\mathbf{R^{(f)}}$ | $M \times M$ | $48 \times 48$ | `cov_mat_freq` |
| $\mathbf{\Pi}_n$ | $M \times K_n$ | $48 \times 24$ | pilot spreading matrix |
| $\mathbf{\Sigma}_n$ | $K_n \times K_n$ | $24 \times 24$ diagonal | LS error at pilots |
| $\mathbf{\Xi}_n$ | $M \times M$ | $48 \times 48$ diagonal | zeros non-pilot subcarrier columns |
| $\mathbf{\Sigma}_n^{(1)}$ | $M \times M$ diag | $48 \times 48$ diag | output error variances |

Symbols without pilots are skipped ($K_n = 0$).

---

### Pass 2 — variance scaling (frequency → time handoff)

After pass 1, rescale each estimate and its variance so the statistics match what the **time** pass expects (uses **R⁽ᶠ⁾**):

**Channel estimate**

$$
\left[\hat{\mathbf{h}}_n^{(2)}\right]_m = s_{n,m} \left[\hat{\mathbf{h}}_n^{(1)}\right]_m
$$

**Error variance**

$$
\left[\mathbf{\Sigma}_n^{(2)}\right]_{m,m} = s_{n,m}\left(s_{n,m}-1\right)\left[\hat{\mathbf{\Sigma}}_n^{(1)}\right]_{m,m} + \left(1-s_{n,m}\right)\left[\mathbf{R^{(f)}}\right]_{m,m} + s_{n,m}\left[\mathbf{\Sigma}_n^{(1)}\right]_{m,m}
$$

**Scaling factor**

$$
s_{n,m} = \frac{2 \left[\mathbf{R^{(f)}}\right]_{m,m}}{\left[\mathbf{R^{(f)}}\right]_{m,m} - \left[\mathbf{\Sigma}_n^{(1)}\right]_{m,m} + \left[\hat{\mathbf{\Sigma}}_n^{(1)}\right]_{m,m}}, \quad \hat{\mathbf{\Sigma}}_n^{(1)} = \mathbf{A}_n \mathbf{R^{(f)}} \mathbf{A}_n^{\mathrm{H}}
$$

Applied on the same DMRS symbols $n \in \{2, 11\}$ as pass 1. Output $\hat{\mathbf{h}}_n^{(2)}$ and $\mathbf{\Sigma}_n^{(2)}$ form the rows of $\hat{\mathbf{H}}^{(2)}$ used as input to pass 3.

| Matrix | Shape | In this script | Meaning |
|--------|-------|----------------|---------|
| $\hat{\mathbf{h}}_n^{(1)}$ (input) | $M \times 1$ | $48 \times 1$ | from pass 1 |
| $\mathbf{\Sigma}_n^{(1)}$ (input) | $M \times M$ diag | $48 \times 48$ diag | from pass 1 |
| $\hat{\mathbf{h}}_n^{(2)}$ (output) | $M \times 1$ | $48 \times 1$ | scaled row for time pass |
| $\mathbf{\Sigma}_n^{(2)}$ (output) | $M \times M$ diag | $48 \times 48$ diag | scaled error variances |
| $\hat{\mathbf{\Sigma}}_n^{(1)}$ | $M \times M$ | $48 \times 48$ | $\mathbf{A}_n \mathbf{R^{(f)}} \mathbf{A}_n^{\mathrm{H}}$ |

---

### Pass 3 — time interpolation (per subcarrier)

For each subcarrier $m = 0,\ldots,M-1$:

**Channel estimate**

$$
\hat{\mathbf{h}}_m^{(3)} = \mathbf{B}_m \tilde{\mathbf{h}}_m^{(2)}
$$

**Error variance** (diagonal; Sionna keeps only the diagonal of the matrix MSE)

$$
\mathbf{\Sigma}_m^{(3)} = \mathrm{diag}\left( \mathbf{R^{(t)}} - \mathbf{B}_m \tilde{\mathbf{\Xi}}_m \mathbf{R^{(t)}} \right)
$$

Equivalently, for each OFDM symbol $n$: $[\mathbf{\Sigma}_m^{(3)}]_{n,n} = [\mathbf{R^{(t)}}]_{n,n} - [\mathbf{B}_m \tilde{\mathbf{\Xi}}_m \mathbf{R^{(t)}}]_{n,n}$.

$\tilde{\mathbf{\Xi}}_m$ is an $N \times N$ diagonal matrix that **zeros columns for OFDM symbols without raw pilots** on subcarrier $m$. Pass 1 fills both DMRS rows for every $m$; pass 2 rescales them. Pass 3 therefore uses $L_m=2$ for all subcarriers (including the 24 comb-gap SC where $L_m^{\mathrm{raw}}=0$).

**Full grid** after pass 3:

$$
\hat{\mathbf{H}}^{(3)} = \left[ \hat{\mathbf{h}}_1^{(3)} \;\cdots\; \hat{\mathbf{h}}_M^{(3)} \right]
$$

**Interpolation matrix**:

$$
\mathbf{B}_m = \bar{\mathbf{B}}_m \tilde{\mathbf{\Pi}}_m^\intercal, \qquad
\bar{\mathbf{B}}_m = \underset{\mathbf{Z} \in \mathbb{C}^{N \times L_m}}{\arg\min} \left\| \mathbf{Z} \left( \tilde{\mathbf{\Pi}}_m^\intercal \mathbf{R^{(t)}} \tilde{\mathbf{\Pi}}_m + \tilde{\mathbf{\Sigma}}_m^{(2)} \right) - \mathbf{R^{(t)}} \tilde{\mathbf{\Pi}}_m \right\|_{\mathrm{F}}^2
$$

| Matrix | Shape | In this script | Meaning |
|--------|-------|----------------|---------|
| $\tilde{\mathbf{h}}_m^{(2)}$ | $N \times 1$ | $14 \times 1$ | column from pass 2 rows |
| $\hat{\mathbf{h}}_m^{(3)}$ (output) | $N \times 1$ | $14 \times 1$ | full column after time interpolation |
| $\bar{\mathbf{B}}_m$ | $N \times L_m$ | $14 \times 2$ | core time-interpolation coefficients |
| $\mathbf{B}_m$ | $N \times N$ | $14 \times 14$ | $\bar{\mathbf{B}}_m \tilde{\mathbf{\Pi}}_m^\intercal$ |
| $\mathbf{R^{(t)}}$ | $N \times N$ | $14 \times 14$ | `cov_mat_time` |
| $\tilde{\mathbf{\Pi}}_m$ | $N \times L_m$ | $14 \times 2$ | $L_m=2$ for all $m$ |
| $\tilde{\mathbf{\Sigma}}_m^{(2)}$ | $L_m \times L_m$ | $2 \times 2$ diagonal | error variances at symbols 2 and 11 |
| $\tilde{\mathbf{\Xi}}_m$ | $N \times N$ | $14 \times 14$ diagonal | zeros non-pilot OFDM-symbol columns |
| $\mathbf{\Sigma}_m^{(3)}$ | $N \times N$ diag | $14 \times 14$ diag | output error variances |

After pass 3, **all** REs have a channel estimate. Our script (`order="f-t"`) stops here; final output is $\hat{\mathbf{H}}^{(3)}$ with `err_var` from $\mathbf{\Sigma}_m^{(3)}$.

---

### Pass 4 — variance scaling (time → spatial handoff; optional)

Only when `order` includes `"s"`. Same pattern as pass 2, but along **time** using **R⁽ᵗ⁾**:

**Channel estimate**

$$
\left[\hat{\mathbf{h}}_m^{(4)}\right]_n = \gamma_{m,n} \left[\hat{\mathbf{h}}_m^{(3)}\right]_n
$$

**Error variance**

$$
\left[\mathbf{\Sigma}_m^{(4)}\right]_{n,n} = \gamma_{m,n}\left(\gamma_{m,n}-1\right)\left[\hat{\mathbf{\Sigma}}_m^{(3)}\right]_{n,n} + \left(1-\gamma_{m,n}\right)\left[\mathbf{R^{(t)}}\right]_{n,n} + \gamma_{m,n}\left[\mathbf{\Sigma}_m^{(3)}\right]_{n,n}
$$

**Scaling factor**

$$
\gamma_{m,n} = \frac{2 \left[\mathbf{R^{(t)}}\right]_{n,n}}{\left[\mathbf{R^{(t)}}\right]_{n,n} - \left[\mathbf{\Sigma}_m^{(3)}\right]_{n,n} + \left[\hat{\mathbf{\Sigma}}_m^{(3)}\right]_{n,n}}, \quad \hat{\mathbf{\Sigma}}_m^{(3)} = \mathbf{B}_m \mathbf{R^{(t)}} \mathbf{B}_m^{\mathrm{H}}
$$

---

### Pass 5 — spatial smoothing (`order` includes `"s"`; optional)

For each resource element $(n,m)$ with $L$ receive antennas (not used in this script):

**Channel estimate**

$$
\hat{\mathbf{h}}^{(5)} = \mathbf{C} \hat{\mathbf{h}}^{(4)}
$$

**Error variance** (diagonal; Sionna keeps only the diagonal of the matrix MSE)

$$
\mathbf{\Sigma}^{(5)} = \mathrm{diag}\left( \mathbf{R^{(s)}} - \mathbf{C}\mathbf{R^{(s)}} \right)
$$

Equivalently, for each antenna $\ell$: $[\mathbf{\Sigma}^{(5)}]_{\ell,\ell} = [\mathbf{R^{(s)}}]_{\ell,\ell} - [\mathbf{C}\mathbf{R^{(s)}}]_{\ell,\ell}$.

**Smoothing matrix** (requires `cov_mat_space`):

$$
\mathbf{C} = \mathbf{R^{(s)}} \left( \mathbf{R^{(s)}} + \mathbf{\Sigma}^{(4)} \right)^{-1}
$$

No scaling is applied after the last pass. Input $\hat{\mathbf{h}}^{(4)}$ and $\mathbf{\Sigma}^{(4)}$ come from pass 4.

| Matrix | Shape | In this script | Meaning |
|--------|-------|----------------|---------|
| $\hat{\mathbf{h}}^{(4)}$ (input) | $L \times 1$ | $1 \times 1$ | channel vector across RX antennas at one RE |
| $\hat{\mathbf{h}}^{(5)}$ (output) | $L \times 1$ | $1 \times 1$ | smoothed channel vector |
| $\mathbf{R^{(s)}}$ | $L \times L$ | $1 \times 1$ | `cov_mat_space` (spatial covariance) |
| $\mathbf{C}$ | $L \times L$ | $1 \times 1$ | spatial smoothing matrix |
| $\mathbf{\Sigma}^{(4)}$ | $L \times L$ | $1 \times 1$ diagonal | error covariance before spatial pass |
| $\mathbf{\Sigma}^{(5)}$ | $L \times L$ | $1 \times 1$ diagonal | error covariance after spatial pass |

Our example uses `order="f-t"` only (passes 1–3; no pass 4 or 5).

---

## Matrix shapes (from `linear_vs_lmmse_interpolator.py`)

Default `PUSCHConfig()` with `dmrs.additional_position = 1` and `num_cdm_groups_without_data = 1`:

| Grid parameter | Symbol | Value in script |
|----------------|--------|-----------------|
| OFDM symbols | $N$ | **14** (`rg.num_ofdm_symbols`) |
| PUSCH allocation | $N_{\mathrm{RB}}$ | **4 RBs** (`pusch_config.num_resource_blocks`) |
| Subcarriers on grid | $M$ | **48** = $N_{\mathrm{RB}} \times 12$ (`rg.fft_size` = `rg.num_effective_subcarriers`; no guard carriers) |
| DMRS symbols (0-based) | — | **2** and **11** (2 DMRS positions) |
| Pilots per DMRS symbol | $K_n$ | **24** = $6 \times N_{\mathrm{RB}}$ on symbols 2, 11; **0** otherwise |
| Raw pilots per subcarrier | $L_m^{\mathrm{raw}}$ | **2** on 24 comb SC; **0** on 24 gap SC (before interpolation) |
| Total pilot REs | $P$ | **48** = $2 \times 24$ (`pilot_pattern.num_pilot_symbols`) |
| RX antennas | $L$ | **1** (SISO) |

### DMRS pilot density

Each PRB has **12 REs** in frequency on one OFDM symbol. With **DMRS config type 1** and **one CDM group**, only **6 REs/PRB** carry DMRS (frequency comb); the other 6 carry PUSCH data.

- $K_n = 6 \times N_{\mathrm{RB}} = 24$ on DMRS symbols $n \in \{2, 11\}$; comb pattern **6/12 REs per PRB** (e.g. SC 0,2,4,… are pilots)
- $K_n < M$ → pass 1 **interpolates in frequency** from 24 pilots to all 48 subcarriers
- $\mathbf{\Pi}_n$ is **$48 \times 24$**; $\bar{\mathbf{A}}_n$ is **$48 \times 24$**

The script sets `pusch_config.dmrs.num_cdm_groups_without_data = 1` explicitly (Sionna default is 2 CDM groups, which uses both frequency combs and doubles DMRS REs per PRB).

### Verified shape summary

| Object | Shape | Notes |
|--------|-------|-------|
| $\hat{\mathbf{H}}$ (one stream) | $14 \times 48$ | $N \times M$ resource grid |
| **R⁽ᶠ⁾** | $48 \times 48$ | `cov_mat_freq` |
| **R⁽ᵗ⁾** | $14 \times 14$ | `cov_mat_time` |
| **Pass 1** (symbols 2, 11 only) | | |
| $\hat{\mathbf{h}}_n$ (input row) | $48 \times 1$ | $K_n=24$ non-zero at pilot SC; zeros on comb gaps |
| $\mathbf{\Pi}_n$ | $48 \times 24$ | |
| $\bar{\mathbf{A}}_n$ | $48 \times 24$ | |
| $\mathbf{A}_n$ | $48 \times 48$ | $\bar{\mathbf{A}}_n \mathbf{\Pi}_n^\intercal$ |
| $\mathbf{\Sigma}_n$ | $24 \times 24$ | diagonal |
| $\mathbf{\Xi}_n$ | $48 \times 48$ | diagonal; zeros non-pilot columns |
| $\mathbf{\Sigma}_n^{(1)}$ | $M \times M$ diag | $\mathrm{diag}(\mathbf{R^{(f)}} - \mathbf{A}_n \mathbf{\Xi}_n \mathbf{R^{(f)}})$ |
| **Pass 2** (variance scaling, symbols 2, 11) | | |
| $\hat{\mathbf{h}}_n^{(2)}$ | $M \times 1$ | scaled row; input to pass 3 |
| $\mathbf{\Sigma}_n^{(2)}$ | $M \times M$ diag | scaled error variances |
| **Pass 3** (time, all $m = 0,\ldots,47$) | | |
| $\tilde{\mathbf{h}}_m^{(2)}$ | $N \times 1$ | column built from pass 2 rows |
| $\tilde{\mathbf{\Pi}}_m$ | $14 \times 2$ | $L_m = 2$ for every $m$ |
| $\bar{\mathbf{B}}_m$ | $14 \times 2$ | |
| $\mathbf{B}_m$ | $14 \times 14$ | |
| $\tilde{\mathbf{\Sigma}}_m^{(2)}$ | $2 \times 2$ diag | error variances at symbols 2 and 11 |
| $\tilde{\mathbf{\Xi}}_m$ | $14 \times 14$ diag | zeros non-pilot OFDM-symbol columns |
| $\mathbf{\Sigma}_m^{(3)}$ | $14 \times 14$ diag | $\mathrm{diag}(\mathbf{R^{(t)}} - \mathbf{B}_m \tilde{\mathbf{\Xi}}_m \mathbf{R^{(t)}})$ |
| **Pass 4** (variance scaling before spatial; optional) | | |
| $\mathbf{\Sigma}_m^{(4)}$ | $14 \times 14$ diag | only if `order` includes `"s"` |
| **Pass 5** (spatial; optional, per RE) | | |
| $\mathbf{\Sigma}^{(4)}$ | $L \times L$ diag | input error before spatial pass |
| $\mathbf{C}$ | $L \times L$ | $\mathbf{R^{(s)}}(\mathbf{R^{(s)}} + \mathbf{\Sigma}^{(4)})^{-1}$ |
| $\mathbf{\Sigma}^{(5)}$ | $L \times L$ diag | $\mathrm{diag}(\mathbf{R^{(s)}} - \mathbf{C}\mathbf{R^{(s)}})$ |

### Sionna tensor shapes

Shapes below omit the **batch** dimension (`batch_size = 16` in the script). The same processing applies independently to each batch item.

| Tensor | Shape | Meaning |
|--------|-------|---------|
| `y` (received grid) | `[1, 1, 14, 48]` | num_rx, num_rx_ant, OFDM symbols, subcarriers |
| `h_hat` (LS at pilots, internal) | `[1, 1, 1, 1, 48]` | num_rx, num_rx_ant, num_tx, num_streams, pilot REs |
| `err_var` (LS error var., internal) | `[1, 1, 1, 1, 48]` | same layout as LS `h_hat` at pilots |
| `h_hat` (output, full grid) | `[1, 1, 1, 1, 14, 48]` | channel estimate on every RE |
| `err_var` (output) | `[1, 1, 1, 1, 14, 48]` | error variance on every RE |
| `h_true` (from channel) | `[1, 1, 1, 1, 14, 48]` | true channel (used for MSE) |
| `cov_mat_freq` | `[48, 48]` | **R⁽ᶠ⁾** |
| `cov_mat_time` | `[14, 14]` | **R⁽ᵗ⁾** |

Mapping to math notation: $\hat{\mathbf{H}} \in \mathbb{C}^{N \times M} = \mathbb{C}^{14 \times 48}$ for one antenna/stream.

---

## MMSE and LMMSE (detailed)

Notation: **x** is the random vector to estimate; **y** is the observation; $\hat{\mathbf{x}}$ is the estimate. Complex vectors use $\|\mathbf{v}\|^2=\mathbf{v}^{\mathrm{H}}\mathbf{v}$.

Four parts: [MMSE](#mmse) → [LMMSE](#lmmse) → LMMSE in [frequency](#lmmse-in-frequency-domain-interpolation) (pass 1), [time](#lmmse-in-time-domain-interpolation) (pass 3), and [spatial](#lmmse-in-spatial-domain-smoothing) (pass 5; optional). Passes 2 and 4 are variance scaling between domains. Full pipeline: [Sionna's algorithm](#sionnas-algorithm-orderf-t).

---

### MMSE

#### Definition

Given random vectors $\mathbf{x}$ (unknown) and $\mathbf{y}$ (observed), the **minimum mean-square error (MMSE)** estimator minimizes average squared error over **all** (possibly nonlinear) functions $f$:

$$
\hat{\mathbf{x}}_{\mathrm{MMSE}} = \underset{f}{\arg\min}\; \mathbb{E}\left\| \mathbf{x} - f(\mathbf{y}) \right\|^2
$$

The minimum is the **Bayes MMSE** (or $J_{\mathrm{MMSE}}$):

$$
J_{\mathrm{MMSE}} = \mathbb{E}\left\| \mathbf{x} - \hat{\mathbf{x}}_{\mathrm{MMSE}} \right\|^2
$$

MMSE is **Bayes estimation** with squared-error loss. No linearity is assumed; the solution depends on the full joint distribution $p(\mathbf{x},\mathbf{y})$.

#### Orthogonality principle

$\hat{\mathbf{x}}=f(\mathbf{y})$ is the MMSE estimator **if and only if** the error is orthogonal to every function of the data:

$$
\mathbb{E}\left[ (\mathbf{x} - f(\mathbf{y}))\, g(\mathbf{y})^{\mathrm{H}} \right] = \mathbf{0}
\qquad \text{for all } g(\cdot)
$$

*Sketch:* perturb $f(\mathbf{y})$ by $\alpha g(\mathbf{y})$; at the optimum, $\frac{\partial}{\partial\alpha}\mathbb{E}\|\mathbf{x}-f(\mathbf{y})-\alpha g(\mathbf{y})\|^2\big|_{\alpha=0}=0$ for every $g$.

#### Conditional mean

The unique MMSE estimator is the **posterior mean**:

$$
\hat{\mathbf{x}}_{\mathrm{MMSE}} = \mathbb{E}[\mathbf{x}\mid\mathbf{y}]
$$

*Sketch:* decompose $\mathbf{x}-f(\mathbf{y})=(\mathbf{x}-\mathbb{E}[\mathbf{x}\mid\mathbf{y}])+(\mathbb{E}[\mathbf{x}\mid\mathbf{y}]-f(\mathbf{y}))$; cross terms vanish when conditioning on $\mathbf{y}$, so any $f\neq\mathbb{E}[\mathbf{x}\mid\mathbf{y}]$ adds excess MSE.

#### MMSE value and error covariance

For any estimator $\hat{\mathbf{x}}$, the **matrix MSE** is

$$
\mathbf{M} = \mathbb{E}\left[ (\mathbf{x}-\hat{\mathbf{x}})(\mathbf{x}-\hat{\mathbf{x}})^{\mathrm{H}} \right]
$$

For the MMSE estimator:

$$
J_{\mathrm{MMSE}} = \mathrm{tr}\,\mathbb{E}\big[\mathrm{Cov}(\mathbf{x}\mid\mathbf{y})\big]
$$

Each diagonal $[\mathbf{M}]_{i,i}$ is the MMSE for estimating $x_i$. **Law of total variance:**

$$
\mathbb{E}\|\mathbf{x}-\hat{\mathbf{x}}\|^2
= \underbrace{\mathbb{E}\|\mathbf{x}-\mathbb{E}[\mathbf{x}\mid\mathbf{y}]\|^2}_{J_{\mathrm{MMSE}}}
+ \underbrace{\mathbb{E}\|\mathbb{E}[\mathbf{x}\mid\mathbf{y}]-\hat{\mathbf{x}}\|^2}_{\text{excess if not MMSE}}
$$

#### When MMSE is linear (jointly $\mathcal{CN}$)

If $\mathbf{x},\mathbf{y}$ are jointly circularly symmetric complex Gaussian with zero mean,

$$
\begin{bmatrix} \mathbf{x} \\ \mathbf{y} \end{bmatrix}
\sim \mathcal{CN}\!\left(\mathbf{0},
\begin{bmatrix} \mathbf{R}_{\mathbf{x}\mathbf{x}} & \mathbf{R}_{\mathbf{x}\mathbf{y}} \\
               \mathbf{R}_{\mathbf{y}\mathbf{x}} & \mathbf{R}_{\mathbf{y}\mathbf{y}} \end{bmatrix}
\right)
$$

then $\mathbb{E}[\mathbf{x}\mid\mathbf{y}]$ is linear in $\mathbf{y}$ (see [LMMSE](#lmmse) below). Wireless channels are modeled as $\mathcal{CN}$, so for channel interpolation the MMSE benchmark reduces to LMMSE.

---

### LMMSE

#### Definition

**Linear MMSE (LMMSE)** restricts the estimator to an **affine** map of the observation (scalar form $\hat{X}_L=aY+b$; vector form $\hat{\mathbf{x}}=\mathbf{A}\mathbf{y}+\mathbf{b}$) and chooses the coefficients to minimize MSE. When $\mathbb{E}[\mathbf{x}]=\mathbb{E}[\mathbf{y}]=\mathbf{0}$ (zero-mean channels in Sionna), the bias term vanishes and $\hat{\mathbf{x}}=\mathbf{A}\mathbf{y}$.

| | MMSE | LMMSE |
|---|------|-------|
| Form | $f(\mathbf{y})$, any $f$ | $\mathbf{A}\mathbf{y}+\mathbf{b}$ only |
| Needs | full $p(\mathbf{x},\mathbf{y})$ | means, $\mathbf{R}_{\mathbf{x}\mathbf{x}}$, $\mathbf{R}_{\mathbf{x}\mathbf{y}}$, $\mathbf{R}_{\mathbf{y}\mathbf{y}}$ |
| Gaussian case | $\mathbb{E}[\mathbf{x}\mid\mathbf{y}]$ | same as MMSE |

#### Proof of the LMMSE estimator

We follow the same method as [Probability Course, §9.1.6](https://www.probabilitycourse.com/chapter9/9_1_6_linear_MMSE_estimat_of_random_vars.php): expand the MSE as a quadratic in the estimator coefficients, set partial derivatives to zero, solve the resulting linear system, then compute the minimum MSE and verify orthogonality.

---

**Theorem (scalar, real).** Let $X,Y$ be real random variables with finite second moments and correlation coefficient $\rho$. Consider

$$
h(a,b)=\mathbb{E}\big[(X-aY-b)^2\big].
$$

Then:

1. $h(a,b)$ is minimized at
   $$
   a^*=\frac{\mathrm{Cov}(X,Y)}{\mathrm{Var}(Y)}, \qquad b^*= \mathbb{E}[X]-a^*\mathbb{E}[Y].
   $$
2. The minimum MSE is
   $$
   h(a^*,b^*)=\big(1-\rho^2\big)\,\mathrm{Var}(X).
   $$
3. The estimation error $\tilde{X}=X-a^*Y-b^*$ satisfies the **orthogonality principle**:
   $$
   \mathbb{E}[\tilde{X}]=0, \qquad \mathbb{E}[\tilde{X}\,Y]=0.
   $$

**Proof.** Expand $h(a,b)$ using linearity of expectation:

$$
\begin{aligned}
h(a,b)
&= \mathbb{E}\big[X^2 + a^2 Y^2 + b^2 - 2aXY - 2bX + 2abY\big] \\
&= \mathbb{E}[X^2] + a^2\mathbb{E}[Y^2] + b^2 - 2a\,\mathbb{E}[XY] - 2b\,\mathbb{E}[X] + 2ab\,\mathbb{E}[Y].
\end{aligned}
$$

Thus $h(a,b)$ is a quadratic in $(a,b)$. Set partial derivatives to zero:

$$
\frac{\partial h}{\partial a} = 2a\,\mathbb{E}[Y^2] + 2b\,\mathbb{E}[Y] - 2\,\mathbb{E}[XY] = 0
\quad\Rightarrow\quad
\mathbb{E}[Y^2]\,a + \mathbb{E}[Y]\,b = \mathbb{E}[XY] \tag{9.4}
$$

$$
\frac{\partial h}{\partial b} = 2b - 2\,\mathbb{E}[X] + 2a\,\mathbb{E}[Y] = 0
\quad\Rightarrow\quad
\mathbb{E}[Y]\,a + b = \mathbb{E}[X] \tag{9.5}
$$

**Solve the $2\times 2$ system.** From (9.5), isolate $b$:

$$
b = \mathbb{E}[X] - \mathbb{E}[Y]\,a. \tag{*}
$$

Substitute into (9.4):

$$
\mathbb{E}[Y^2]\,a + \mathbb{E}[Y]\big(\mathbb{E}[X] - \mathbb{E}[Y]\,a\big) = \mathbb{E}[XY]
$$

Collect terms in $a$:

$$
\big(\mathbb{E}[Y^2] - \mathbb{E}[Y]^2\big)\,a = \mathbb{E}[XY] - \mathbb{E}[X]\,\mathbb{E}[Y].
$$

Recognize variance and covariance:

$$
\mathrm{Var}(Y)\,a = \mathrm{Cov}(X,Y)
\quad\Rightarrow\quad
\boxed{a^* = \frac{\mathrm{Cov}(X,Y)}{\mathrm{Var}(Y)}}.
$$

Plug $a^*$ back into $(*)$:

$$
\boxed{b^* = \mathbb{E}[X] - a^*\,\mathbb{E}[Y]
= \mathbb{E}[X] - \frac{\mathrm{Cov}(X,Y)}{\mathrm{Var}(Y)}\,\mathbb{E}[Y]}.
$$

The Hessian $\mathbf{H}=\begin{bmatrix}2\mathbb{E}[Y^2] & 2\mathbb{E}[Y] \\ 2\mathbb{E}[Y] & 2\end{bmatrix}$ is positive semidefinite (positive definite when $\mathrm{Var}(Y)>0$), so this critical point is the global minimum.

**Minimum MSE.** Equation (9.5) implies $\mathbb{E}[X-a^*Y-b^*]=0$, so at the optimum the MSE equals the error variance:

$$
\begin{aligned}
h(a^*,b^*)
&= \mathrm{Var}(X-a^*Y-b^*) = \mathrm{Var}(X-a^*Y) \\
&= \mathrm{Var}(X) + {a^*}^2\mathrm{Var}(Y) - 2a^*\mathrm{Cov}(X,Y) \\
&= \mathrm{Var}(X) - \frac{\mathrm{Cov}(X,Y)^2}{\mathrm{Var}(Y)} = \big(1-\rho^2\big)\,\mathrm{Var}(X).
\end{aligned}
$$

**Orthogonality.** Substitute the optimal $(a^*,b^*)$ into $\mathbb{E}[(X-a^*Y-b^*)Y]$:

$$
\mathbb{E}[(X-a^*Y-b^*)Y] = \mathbb{E}[XY] - a^*\mathbb{E}[Y^2] - b^*\mathbb{E}[Y] = 0 \quad \text{(by (9.4))}.
$$

Equivalently, $\mathrm{Cov}(\tilde{X},Y)=\mathbb{E}[\tilde{X}\,Y]=0$.

**Estimator form.** The linear MMSE estimate can be written as

$$
\hat{X}_L = a^*(Y-\mathbb{E}[Y]) + \mathbb{E}[X]
= \frac{\mathrm{Cov}(X,Y)}{\mathrm{Var}(Y)}\,(Y-\mathbb{E}[Y]) + \mathbb{E}[X].
$$

**Numerical example.** Let $X \sim \mathrm{Uniform}(1,2)$ and, given $X=x$, let $Y \sim \mathrm{Exponential}(\lambda=1/x)$. We **know the model** (prior on $X$, law of $Y\mid X$) and use it to build $\hat{X}_L$. We do **not** know the realized $X$ for a given trial — only an observation $Y=y$.

**What is $Y$? Where is the exponential used?**

Each trial has two stages:

1. Nature draws the hidden $X$ (uniform on $[1,2]$).
2. Given that value $X=x$, nature draws one **observation** $Y$ from an exponential distribution.

So **$Y$ is what we measure**; **$X$ is what we infer**. The exponential appears only in step 2 — it models how noisy the observation is **given** the true $x$.

**What does $\lambda = 1/x$ mean?**

Write $Y \sim \mathrm{Exponential}(\lambda)$ with **rate** $\lambda > 0$ (same convention as [Probability Course §9.1.6](https://www.probabilitycourse.com/chapter9/9_1_6_linear_MMSE_estimat_of_random_vars.php)). Then

| Quantity | Formula |
|----------|---------|
| PDF | $f_Y(y)=\lambda e^{-\lambda y}$ for $y \ge 0$ |
| Mean | $\mathbb{E}[Y]=1/\lambda$ |
| Variance | $\mathrm{Var}(Y)=1/\lambda^2$ |
| Second moment | $\mathbb{E}[Y^2]=2/\lambda^2$ |

In our example, $\lambda = 1/x$ **depends on the hidden** $x$:

$$
Y \mid X=x \sim \mathrm{Exponential}(\lambda = 1/x)
\quad\Rightarrow\quad
f_{Y\mid X}(y\mid x) = \frac{1}{x}\, e^{-y/x}, \quad y \ge 0.
$$

Substitute $\lambda = 1/x$:

$$
\mathbb{E}[Y \mid X=x] = \frac{1}{\lambda} = x, \qquad
\mathrm{Var}(Y \mid X=x) = \frac{1}{\lambda^2} = x^2, \qquad
\mathbb{E}[Y^2 \mid X=x] = \frac{2}{\lambda^2} = 2x^2.
$$

So **on average, the observation equals the hidden value** ($\mathbb{E}[Y\mid X=x]=x$), but each $Y$ is random and can be much larger or smaller than $x$. That is why we need estimation instead of setting $\hat{X}_L = Y$.

*Concrete draw.* If $X=1.5$, then $\lambda = 1/1.5 = 2/3$ and $Y$ is exponential with mean $1.5$. One possible outcome is $Y=2.0$ (used in Step 4 below) — plausible for $X=1.5$, but not equal to $X$. LMMSE maps $Y=2.0 \mapsto \hat{X}_L \approx 1.52$.

*Step 1 — first and second moments.*

**Why $\mathbb{E}[X]=\frac{3}{2}$?** For $X \sim \mathrm{Uniform}(a,b)$, all values in $[a,b]$ are equally likely, so the mean is the midpoint:

$$
\mathbb{E}[X]=\frac{a+b}{2}=\frac{1+2}{2}=\frac{3}{2}.
$$

Equivalently, with PDF $f_X(x)=1$ on $[1,2]$:

$$
\mathbb{E}[X]=\int_1^2 x\,dx=\left[\frac{x^2}{2}\right]_1^2=\frac{4}{2}-\frac{1}{2}=\frac{3}{2}.
$$

This is the **prior mean** — our guess for $X$ before seeing $Y$. It enters $b^*$ and pulls $\hat{X}_L$ toward $1.5$ when $Y$ is noisy.

**Other moments** (law of iterated expectations — average over both $X$ and $Y$):

$$
\mathbb{E}[Y]=\mathbb{E}[\mathbb{E}[Y\mid X]]=\mathbb{E}[X]=\frac{3}{2}.
$$

The inner step uses $\mathbb{E}[Y\mid X=x]=x$ from the exponential with $\lambda=1/x$.

$$
\mathbb{E}[Y^2]=\mathbb{E}[\mathbb{E}[Y^2\mid X]]=\mathbb{E}[2X^2]=\int_1^2 2x^2\,dx=\frac{14}{3}.
$$

Here $\mathbb{E}[Y^2\mid X=x]=2x^2$ because $\mathbb{E}[Y^2]=2/\lambda^2=2x^2$ when $\lambda=1/x$.

$$
\mathbb{E}[XY]=\mathbb{E}[X\,\mathbb{E}[Y\mid X]]=\mathbb{E}[X^2]=\int_1^2 x^2\,dx=\frac{7}{3}.
$$

*Step 2 — variance and covariance.*

For $\mathrm{Uniform}(1,2)$, $\mathrm{Var}(X)=\frac{(b-a)^2}{12}=\frac{1}{12}$. Then:

$$
\mathrm{Var}(Y)=\mathbb{E}[Y^2]-(\mathbb{E}[Y])^2=\frac{14}{3}-\frac{9}{4}=\frac{29}{12}, \qquad
\mathrm{Var}(X)=\frac{1}{12}.
$$

$$
\mathrm{Cov}(X,Y)=\mathbb{E}[XY]-\mathbb{E}[X]\mathbb{E}[Y]=\frac{7}{3}-\frac{9}{4}=\frac{1}{12}.
$$

*Step 3 — LMMSE coefficients (theorem part 1).*

$$
a^*=\frac{\mathrm{Cov}(X,Y)}{\mathrm{Var}(Y)}
=\frac{1/12}{29/12}=\frac{1}{29}, \qquad
b^*=\mathbb{E}[X]-a^*\mathbb{E}[Y]=\frac{3}{2}-\frac{1}{29}\cdot\frac{3}{2}=\frac{42}{29}.
$$

So

$$
\hat{X}_L=\frac{1}{29}\,Y+\frac{42}{29}
=\frac{1}{29}(Y-\mathbb{E}[Y])+\mathbb{E}[X].
$$

*Step 4 — one observed value.* Suppose one trial yields **$Y=2.0$** (we do not see $X$). Plug in:

$$
\hat{X}_L = \frac{1}{29}(2.0) + \frac{42}{29} = \frac{44}{29} \approx 1.52.
$$

| Quantity | Value |
|----------|-------|
| Prior mean $\mathbb{E}[X]$ | $1.50$ |
| Observation $Y$ | $2.00$ |
| LMMSE estimate $\hat{X}_L$ | $\approx 1.52$ |

The observation $Y=2.0$ is above the prior mean $1.5$, so $\hat{X}_L > 1.5$, but only slightly — $a^*=1/29$ **shrinks** the estimate toward $\mathbb{E}[X]$ because $Y$ is a noisy proxy for $X$ ($\mathrm{Var}(Y)=29/12 \gg \mathrm{Var}(X)=1/12$). We do **not** set $\hat{X}_L = Y$.

*Step 5 — MSE (theorem part 2).*

$$
\rho^2=\frac{\mathrm{Cov}(X,Y)^2}{\mathrm{Var}(X)\,\mathrm{Var}(Y)}
=\frac{(1/12)^2}{(1/12)(29/12)}=\frac{1}{29}.
$$

$$
h(a^*,b^*)=(1-\rho^2)\,\mathrm{Var}(X)
=\left(1-\frac{1}{29}\right)\frac{1}{12}=\frac{7}{87}.
$$

*Step 6 — orthogonality check (theorem part 3).* With $\tilde{X}=X-\hat{X}_L=X-\frac{Y}{29}-\frac{42}{29}$:

$$
\mathbb{E}[\tilde{X}\,Y]
=\mathbb{E}[XY]-\frac{\mathbb{E}[Y^2]}{29}-\frac{42}{29}\mathbb{E}[Y]
=\frac{7}{3}-\frac{14}{3\cdot 29}-\frac{42}{29}\cdot\frac{3}{2}=0.
$$

Also $\mathbb{E}[\tilde{X}]=\mathbb{E}[X]-a^*\mathbb{E}[Y]-b^*=0$ by construction of $b^*$ from (9.5).

---

**Theorem (vector, complex).** Let $\mathbf{x},\mathbf{y}\in\mathbb{C}^n$ have finite second moments. Consider the affine estimator $\hat{\mathbf{x}}=\mathbf{A}\mathbf{y}+\mathbf{b}$ and

$$
h(\mathbf{A},\mathbf{b})=\mathbb{E}\big[(\mathbf{x}-\mathbf{A}\mathbf{y}-\mathbf{b})^{\mathrm{H}}(\mathbf{x}-\mathbf{A}\mathbf{y}-\mathbf{b})\big].
$$

Then:

1. $h$ is minimized at
   $$
   \mathbf{A}^*\mathbf{R}_{\mathbf{y}\mathbf{y}} = \mathbf{R}_{\mathbf{x}\mathbf{y}}, \qquad
   \mathbf{b}^* = \mathbb{E}[\mathbf{x}] - \mathbf{A}^*\mathbb{E}[\mathbf{y}],
   $$
   and when $\mathbf{R}_{\mathbf{y}\mathbf{y}}$ is invertible, $\mathbf{A}^*=\mathbf{R}_{\mathbf{x}\mathbf{y}}\mathbf{R}_{\mathbf{y}\mathbf{y}}^{-1}$.
2. The error covariance at the optimum is
   $$
   \mathbf{M}=\mathbb{E}[(\mathbf{x}-\mathbf{A}^*\mathbf{y}-\mathbf{b}^*)(\cdot)^{\mathrm{H}}]
   = \mathbf{R}_{\mathbf{x}\mathbf{x}} - \mathbf{R}_{\mathbf{x}\mathbf{y}}\mathbf{R}_{\mathbf{y}\mathbf{y}}^{-1}\mathbf{R}_{\mathbf{y}\mathbf{x}}
   $$
   (for zero means, this is the matrix analogue of $(1-\rho^2)\mathrm{Var}(X)$).
3. The error $\tilde{\mathbf{x}}=\mathbf{x}-\mathbf{A}^*\mathbf{y}-\mathbf{b}^*$ satisfies
   $$
   \mathbb{E}[\tilde{\mathbf{x}}]=\mathbf{0}, \qquad \mathbb{E}[\tilde{\mathbf{x}}\mathbf{y}^{\mathrm{H}}]=\mathbf{0}.
   $$

**Proof (same steps as scalar).** Expand $\mathbf{e}=\mathbf{x}-\mathbf{A}\mathbf{y}-\mathbf{b}$:

$$
h(\mathbf{A},\mathbf{b})
= \mathbb{E}[\mathbf{e}^{\mathrm{H}}\mathbf{e}]
= \mathbb{E}[\mathbf{x}^{\mathrm{H}}\mathbf{x}]
  - 2\,\mathrm{Re}\,\mathbb{E}[\mathbf{x}^{\mathrm{H}}\mathbf{A}\mathbf{y}]
  - 2\,\mathrm{Re}\,\mathbb{E}[\mathbf{b}^{\mathrm{H}}\mathbf{x}]
  + \mathbb{E}[\mathbf{y}^{\mathrm{H}}\mathbf{A}^{\mathrm{H}}\mathbf{A}\mathbf{y}]
  + 2\,\mathrm{Re}\,\mathbb{E}[\mathbf{b}^{\mathrm{H}}\mathbf{A}\mathbf{y}]
  + \mathbf{b}^{\mathrm{H}}\mathbf{b}.
$$

This is a quadratic in $\mathbf{A}$ and $\mathbf{b}$. Setting gradients to zero (w.r.t. $\mathbf{b}^*$ and $\mathbf{A}^*$) gives the vector analogues of (9.5) and (9.4):

$$
\mathbf{b} = \mathbb{E}[\mathbf{x}] - \mathbf{A}\,\mathbb{E}[\mathbf{y}], \tag{9.5'}
$$

$$
\mathbf{A}\,\mathbf{R}_{\mathbf{y}\mathbf{y}} = \mathbf{R}_{\mathbf{x}\mathbf{y}}. \tag{9.4'}
$$

Solve $(9.4')$ for $\mathbf{A}^*$, then $(9.5')$ for $\mathbf{b}^*$.

**Minimum MSE / error covariance.** At the optimum, $(9.5')$ gives $\mathbb{E}[\tilde{\mathbf{x}}]=\mathbf{0}$, so

$$
h(\mathbf{A}^*,\mathbf{b}^*) = \mathbb{E}[\tilde{\mathbf{x}}^{\mathrm{H}}\tilde{\mathbf{x}}] = \mathrm{tr}(\mathbf{M}).
$$

Expanding $\mathbb{E}[\tilde{\mathbf{x}}\tilde{\mathbf{x}}^{\mathrm{H}}]$ and using $\mathbf{A}^*\mathbf{R}_{\mathbf{y}\mathbf{y}}=\mathbf{R}_{\mathbf{x}\mathbf{y}}$ yields the Schur-complement form of $\mathbf{M}$ above — the same algebra as the scalar step $\mathrm{Var}(X)-\mathrm{Cov}(X,Y)^2/\mathrm{Var}(Y)$.

**Orthogonality.** From $(9.4')$:

$$
\mathbb{E}[(\mathbf{x}-\mathbf{A}^*\mathbf{y}-\mathbf{b}^*)\mathbf{y}^{\mathrm{H}}]=\mathbf{0}.
$$

**Zero-mean case (Sionna).** When $\mathbb{E}[\mathbf{x}]=\mathbb{E}[\mathbf{y}]=\mathbf{0}$, $(9.5')$ gives $\mathbf{b}^*=\mathbf{0}$ and

$$
\boxed{\hat{\mathbf{x}}_{\mathrm{LMMSE}} = \mathbf{R}_{\mathbf{x}\mathbf{y}}\mathbf{R}_{\mathbf{y}\mathbf{y}}^{-1}\mathbf{y}}.
$$

For jointly Gaussian $\mathbf{x},\mathbf{y}$, this equals $\mathbb{E}[\mathbf{x}\mid\mathbf{y}]$ — LMMSE and MMSE coincide.

**Summary**

| Quantity | Scalar | Vector (zero mean) |
|----------|--------|---------------------|
| Estimator | $\frac{\mathrm{Cov}(X,Y)}{\mathrm{Var}(Y)}(Y-\mathbb{E}Y)+\mathbb{E}X$ | $\mathbf{R}_{\mathbf{x}\mathbf{y}}\mathbf{R}_{\mathbf{y}\mathbf{y}}^{-1}\mathbf{y}$ |
| Normal equations | (9.4), (9.5) | $\mathbf{A}\mathbf{R}_{\mathbf{y}\mathbf{y}}=\mathbf{R}_{\mathbf{x}\mathbf{y}}$ |
| Orthogonality | $\mathbb{E}[\tilde{X}Y]=0$ | $\mathbb{E}[(\mathbf{x}-\mathbf{A}\mathbf{y})\mathbf{y}^{\mathrm{H}}]=\mathbf{0}$ |
| MSE | $(1-\rho^2)\mathrm{Var}(X)$ | $\mathbf{R}_{\mathbf{x}\mathbf{x}}-\mathbf{R}_{\mathbf{x}\mathbf{y}}\mathbf{R}_{\mathbf{y}\mathbf{y}}^{-1}\mathbf{R}_{\mathbf{y}\mathbf{x}}$ |

#### Linear measurement model

Often $\mathbf{y}=\mathbf{G}\mathbf{x}+\mathbf{n}$ with $\mathbf{n}\sim\mathcal{CN}(\mathbf{0},\mathbf{R}_{\mathbf{n}\mathbf{n}})$ independent of $\mathbf{x}$. Then

$$
\mathbf{R}_{\mathbf{x}\mathbf{y}} = \mathbf{R}_{\mathbf{x}\mathbf{x}}\mathbf{G}^{\mathrm{H}}, \qquad
\mathbf{R}_{\mathbf{y}\mathbf{y}} = \mathbf{G}\mathbf{R}_{\mathbf{x}\mathbf{x}}\mathbf{G}^{\mathrm{H}} + \mathbf{R}_{\mathbf{n}\mathbf{n}}
$$

and

$$
\mathbf{A} = \mathbf{R}_{\mathbf{x}\mathbf{x}}\mathbf{G}^{\mathrm{H}}
\big(\mathbf{G}\mathbf{R}_{\mathbf{x}\mathbf{x}}\mathbf{G}^{\mathrm{H}} + \mathbf{R}_{\mathbf{n}\mathbf{n}}\big)^{-1}
$$

When observations lie on a **sparse subset** of the full index set, $\mathbf{G}=\mathbf{\Pi}^{\intercal}$ picks/spreads $K$ pilot values onto a length-$D$ grid (see frequency interpolation below). This is the Wiener filter — Step 3 above with $\mathbf{R}_{\mathbf{x}\mathbf{y}}=\mathbf{R}_{\mathbf{x}\mathbf{x}}\mathbf{G}^{\mathrm{H}}$ and $\mathbf{R}_{\mathbf{y}\mathbf{y}}=\mathbf{G}\mathbf{R}_{\mathbf{x}\mathbf{x}}\mathbf{G}^{\mathrm{H}}+\mathbf{R}_{\mathbf{n}\mathbf{n}}$.

#### LMMSE error covariance

From the vector theorem above (minimum MSE step):

$$
\mathbf{M} = \mathbb{E}\left[ (\mathbf{x}-\mathbf{A}\mathbf{y})(\mathbf{x}-\mathbf{A}\mathbf{y})^{\mathrm{H}} \right]
= \mathbf{R}_{\mathbf{x}\mathbf{x}} - \mathbf{R}_{\mathbf{x}\mathbf{y}}\,\mathbf{R}_{\mathbf{y}\mathbf{y}}^{-1}\,\mathbf{R}_{\mathbf{y}\mathbf{x}}
$$

For jointly Gaussian $\mathbf{x},\mathbf{y}$, this equals $\mathrm{Cov}(\mathbf{x}\mid\mathbf{y})$.

#### Scalar example

Estimate $x_1$ from one noisy observation $y_0$ (zero mean). With $R_{ij}=[\mathbf{R}_{\mathbf{x}\mathbf{x}}]_{i,j}$, observation variance $\mathrm{Var}(Y_0)=R_{00}+\sigma^2$, and $\mathrm{Cov}(X_1,Y_0)=R_{01}$, the scalar formula gives

$$
\hat{x}_1 = \frac{\mathrm{Cov}(X_1,Y_0)}{\mathrm{Var}(Y_0)}\,y_0
= \frac{R_{01}}{R_{00} + \sigma^2}\,y_0
$$

Strong correlation and low noise → weight near $R_{01}/R_{00}$; very noisy observation → weight $\to 0$ (shrinks to prior mean $0$). The MSE is $(1-\rho^2)\mathrm{Var}(X_1)$ with $\rho^2=R_{01}^2/\big(R_{00}(R_{00}+\sigma^2)\big)$.

#### Numerical solve (COD)

When $\mathbf{R}_{\mathbf{y}\mathbf{y}}$ is ill-conditioned, Sionna avoids a direct inverse and finds $\mathbf{A}$ via **complete orthogonal decomposition (COD)** — minimum-norm least squares on the normal equations. Same LMMSE solution when the system is well posed.

---

### LMMSE in frequency-domain interpolation

Sionna **pass 1** applies LMMSE along **frequency** on each DMRS OFDM symbol: estimate the full channel row $\mathbf{x}\in\mathbb{C}^{M}$ from $K$ noisy pilot LS values $\mathbf{y}\in\mathbb{C}^{K}$. See [Pass 1](#pass-1--frequency-interpolation-per-ofdm-symbol) for Sionna matrix names ($\mathbf{A}_n$, $\mathbf{\Pi}_n$, …).

#### Problem setup (one OFDM symbol $n$)

| General | Frequency pass (Sionna) | In our script |
|---------|-------------------------|---------------|
| $\mathbf{x}\in\mathbb{C}^{M}$ | channel on $M$ subcarriers ($\mathbf{h}_n$) | $M=48$ |
| $\mathbf{y}\in\mathbb{C}^{K}$ | pilot LS values | $K_n=24$ on symbols 2, 11 |
| $\mathbf{R}_{\mathbf{x}\mathbf{x}}$ | **R⁽ᶠ⁾** (`cov_mat_freq`) | $48\times 48$ |
| $\mathbf{R}_{\mathbf{n}\mathbf{n}}$ | **Σ**$_n$ (diagonal, from `err_var`) | $24\times 24$ |
| $\mathbf{G}=\mathbf{\Pi}_n^{\intercal}$ | pilot spreading / selection | $\mathbf{\Pi}_n$ is $48\times 24$ |

**Measurement model:**

$$
\mathbf{y} = \mathbf{\Pi}_n^{\intercal}\mathbf{x} + \mathbf{n}, \qquad
\mathbf{n} \sim \mathcal{CN}(\mathbf{0}, \mathbf{\Sigma}_n)
$$

LS places the $K$ values on the grid (zeros at non-pilots): $\mathbf{y}_{\mathrm{grid}}=\mathbf{\Pi}_n\mathbf{y}$. Only comb subcarriers carry DMRS (6 REs/PRB with 1 CDM group); pass 1 **interpolates** to all $M$ subcarriers.

#### LMMSE estimate (Wiener filter on pilots)

Substitute $\mathbf{G}=\mathbf{\Pi}_n^{\intercal}$, $\mathbf{R}_{\mathbf{x}\mathbf{x}}=\mathbf{R^{(f)}}$ into the linear measurement model:

$$
\hat{\mathbf{x}}_{\mathrm{LMMSE}}
= \mathbf{R^{(f)}}\mathbf{\Pi}_n
\big(\mathbf{\Pi}_n^{\intercal}\mathbf{R^{(f)}}\mathbf{\Pi}_n + \mathbf{\Sigma}_n\big)^{-1}\mathbf{y}
$$

Define core coefficients $\bar{\mathbf{A}}_n$ ($M\times K_n$) and full-row map $\mathbf{A}_n$ ($M\times M$):

$$
\bar{\mathbf{A}}_n = \mathbf{R^{(f)}}\mathbf{\Pi}_n
\big(\mathbf{\Pi}_n^{\intercal}\mathbf{R^{(f)}}\mathbf{\Pi}_n + \mathbf{\Sigma}_n\big)^{-1}, \qquad
\mathbf{A}_n = \bar{\mathbf{A}}_n\mathbf{\Pi}_n^{\intercal}
$$

Then (matching Sionna pass 1):

$$
\hat{\mathbf{x}} = \bar{\mathbf{A}}_n\mathbf{y} = \mathbf{A}_n\mathbf{y}_{\mathrm{grid}}
\qquad\text{(Sionna: }\hat{\mathbf{h}}_n^{(1)} = \mathbf{A}_n\hat{\mathbf{h}}_n\text{)}
$$

| Matrix | Role |
|--------|------|
| $\mathbf{R^{(f)}}\mathbf{\Pi}_n$ | cross-covariance $\mathbf{R}_{\mathbf{x}\mathbf{y}}$ |
| $\mathbf{\Pi}_n^{\intercal}\mathbf{R^{(f)}}\mathbf{\Pi}_n + \mathbf{\Sigma}_n$ | observation covariance $\mathbf{R}_{\mathbf{y}\mathbf{y}}$ |
| $\bar{\mathbf{A}}_n$ | Wiener weights on pilot vector **y** |
| $\mathbf{A}_n$ | map sparse LS row → full channel row **x** |

Sionna computes $\bar{\mathbf{A}}_n$ by COD:

$$
\bar{\mathbf{A}}_n = \underset{\mathbf{Z}}{\arg\min}
\left\|\mathbf{Z}\big(\mathbf{\Pi}_n^{\intercal}\mathbf{R^{(f)}}\mathbf{\Pi}_n + \mathbf{\Sigma}_n\big)
- \mathbf{R^{(f)}}\mathbf{\Pi}_n\right\|_{\mathrm{F}}^2
$$

#### Error variance after frequency LMMSE

Matrix MSE (before Sionna's diagonal approximation):

$$
\mathbf{M}_n = \mathbf{R^{(f)}} - \mathbf{A}_n\,\mathbf{\Xi}_n\,\mathbf{R^{(f)}}
$$

**Ξ**$_n$ is $M\times M$ diagonal and zeros columns for subcarriers **without raw pilots** on symbol $n$. Output per-subcarrier variances:

$$
\mathbf{\Sigma}_n^{(1)} = \mathrm{diag}(\mathbf{M}_n)
$$

These feed pass 2 (variance scaling) and eventually the time pass. Sionna keeps only $\mathrm{diag}(\mathbf{M}_n)$ — an approximate LMMSE when errors are correlated across subcarriers (see [assumptions](#important-assumptions-from-sionna-docs)).

#### Intuition vs linear interpolation

At each non-pilot subcarrier $m$, $\hat{x}_m$ is a weighted sum of the $K$ pilot LS values. Weights come from **R⁽ᶠ⁾** (how correlated subcarriers are) and **Σ**$_n$ (how noisy each pilot is). **Linear** interpolation uses fixed weights and ignores both — the gap shown in `mse_comparison.png`, especially at large delay spread where **R⁽ᶠ⁾** structure matters most.

---

### LMMSE in time-domain interpolation

Sionna **pass 3** applies LMMSE along **time** on each subcarrier: estimate the full column $\mathbf{x}\in\mathbb{C}^{N}$ (channel across OFDM symbols) from $L_m$ noisy values $\mathbf{y}\in\mathbb{C}^{L_m}$ at DMRS symbol positions. See [Pass 3](#pass-3--time-interpolation-per-subcarrier) for Sionna names ($\mathbf{B}_m$, $\tilde{\mathbf{\Pi}}_m$, …).

**Prerequisite:** pass 1 fills every subcarrier row; pass 2 rescales rows/variances for the time pass. After pass 2, every $m$ has estimates at both DMRS symbols, so $L_m=2$ for all subcarriers (including comb-gap SC where raw pilots were zero before pass 1).

#### Problem setup (one subcarrier $m$)

| General | Time pass (Sionna) | In our script |
|---------|-------------------|---------------|
| $\mathbf{x}\in\mathbb{C}^{N}$ | channel on $N$ OFDM symbols ($\mathbf{h}_m$) | $N=14$ |
| $\mathbf{y}\in\mathbb{C}^{L_m}$ | values at DMRS symbols (from pass 2) | $L_m=2$ (symbols 2, 11) |
| $\mathbf{R}_{\mathbf{x}\mathbf{x}}$ | **R⁽ᵗ⁾** (`cov_mat_time`) | $14\times 14$ |
| $\mathbf{R}_{\mathbf{n}\mathbf{n}}$ | $\tilde{\mathbf{\Sigma}}_m^{(2)}$ (diagonal) | $2\times 2$ |
| $\mathbf{G}=\tilde{\mathbf{\Pi}}_m^{\intercal}$ | pilot spreading / selection | $\tilde{\mathbf{\Pi}}_m$ is $14\times 2$ |

**Measurement model:**

$$
\mathbf{y} = \tilde{\mathbf{\Pi}}_m^{\intercal}\mathbf{x} + \mathbf{n}, \qquad
\mathbf{n} \sim \mathcal{CN}(\mathbf{0}, \tilde{\mathbf{\Sigma}}_m^{(2)})
$$

Input on the symbol grid: $\mathbf{y}_{\mathrm{grid}}=\tilde{\mathbf{\Pi}}_m\mathbf{y}$ (non-DMRS symbols zero before interpolation).

#### LMMSE estimate (Wiener filter)

Same structure as the frequency pass, with **R⁽ᵗ⁾** replacing **R⁽ᶠ⁾**:

$$
\hat{\mathbf{x}}_{\mathrm{LMMSE}}
= \mathbf{R^{(t)}}\tilde{\mathbf{\Pi}}_m
\big(\tilde{\mathbf{\Pi}}_m^{\intercal}\mathbf{R^{(t)}}\tilde{\mathbf{\Pi}}_m + \tilde{\mathbf{\Sigma}}_m^{(2)}\big)^{-1}\mathbf{y}
$$

$$
\bar{\mathbf{B}}_m = \mathbf{R^{(t)}}\tilde{\mathbf{\Pi}}_m
\big(\tilde{\mathbf{\Pi}}_m^{\intercal}\mathbf{R^{(t)}}\tilde{\mathbf{\Pi}}_m + \tilde{\mathbf{\Sigma}}_m^{(2)}\big)^{-1}, \qquad
\mathbf{B}_m = \bar{\mathbf{B}}_m\tilde{\mathbf{\Pi}}_m^{\intercal}
$$

Then:

$$
\hat{\mathbf{x}} = \bar{\mathbf{B}}_m\mathbf{y} = \mathbf{B}_m\mathbf{y}_{\mathrm{grid}}
\qquad\text{(Sionna: }\hat{\mathbf{h}}_m^{(3)} = \mathbf{B}_m\tilde{\mathbf{h}}_m^{(2)}\text{)}
$$

| Matrix | Role |
|--------|------|
| $\mathbf{R^{(t)}}\tilde{\mathbf{\Pi}}_m$ | cross-covariance $\mathbf{R}_{\mathbf{x}\mathbf{y}}$ |
| $\tilde{\mathbf{\Pi}}_m^{\intercal}\mathbf{R^{(t)}}\tilde{\mathbf{\Pi}}_m + \tilde{\mathbf{\Sigma}}_m^{(2)}$ | observation covariance $\mathbf{R}_{\mathbf{y}\mathbf{y}}$ |
| $\bar{\mathbf{B}}_m$ | Wiener weights on **y** |
| $\mathbf{B}_m$ | map sparse symbol column → full column **x** |

COD minimization (same form as pass 1, with **R⁽ᵗ⁾** and $\tilde{\mathbf{\Pi}}_m$).

#### Error variance after time LMMSE

$$
\mathbf{M}_m = \mathbf{R^{(t)}} - \mathbf{B}_m\,\tilde{\mathbf{\Xi}}_m\,\mathbf{R^{(t)}}, \qquad
\mathbf{\Sigma}_m^{(3)} = \mathrm{diag}(\mathbf{M}_m)
$$

**Ξ̃**$_m$ zeros columns for OFDM symbols **without raw pilots** on subcarrier $m$ (symbols other than 2 and 11 before interpolation). Our script (`order="f-t"`) stops here; final grid is $\hat{\mathbf{H}}^{(3)}$ with `err_var` from $\mathbf{\Sigma}_m^{(3)}$.

#### Intuition vs linear interpolation

Between DMRS symbols 2 and 11, **linear** interpolation draws a straight line in time. **LMMSE** uses **R⁽ᵗ⁾**, which shrinks as UE speed rises (Doppler → faster channel variation → weaker time correlation). At high speed, LMMSE trusts pilots more locally; linear over-smooths across the 9-symbol gap.

---

### LMMSE in spatial-domain smoothing

Sionna **pass 5** applies LMMSE across **receive antennas** at each resource element $(n,m)$ when `order` includes `"s"`. Not used in our script (SISO). See [Pass 5](#pass-5--spatial-smoothing-order-includes-s).

Unlike frequency/time passes, **every antenna carries an input estimate** after pass 4 — there is no sparse **Π**. The problem is **denoising / smoothing** a vector $\mathbf{x}\in\mathbb{C}^{L}$ (channel across $L$ RX antennas) from a noisy observation $\mathbf{y}\in\mathbb{C}^{L}$ (pass-4 estimate $\hat{\mathbf{h}}^{(4)}$).

#### Problem setup (one RE)

| General | Spatial pass (Sionna) | Typical MIMO |
|---------|----------------------|--------------|
| $\mathbf{x}\in\mathbb{C}^{L}$ | true channel across RX antennas | $L$ = num_rx_ant |
| $\mathbf{y}\in\mathbb{C}^{L}$ | estimate from pass 4 | same dimension |
| $\mathbf{R}_{\mathbf{x}\mathbf{x}}$ | **R⁽ˢ⁾** (`cov_mat_space`) | $L\times L$ |
| $\mathbf{R}_{\mathbf{n}\mathbf{n}}$ | **Σ**$^{(4)}$ (diagonal, from pass 4) | $L\times L$ diag |
| $\mathbf{G}=\mathbf{I}_L$ | full observation (no sparsity) | — |

**Measurement model:**

$$
\mathbf{y} = \mathbf{x} + \mathbf{n}, \qquad
\mathbf{n} \sim \mathcal{CN}(\mathbf{0}, \mathbf{\Sigma}^{(4)})
$$

Pass 4 scaling aligns $\mathbf{y}$ statistics with **R⁽ˢ⁾** before smoothing.

#### LMMSE estimate (spatial Wiener filter)

With $\mathbf{G}=\mathbf{I}$:

$$
\hat{\mathbf{x}}_{\mathrm{LMMSE}} = \mathbf{R^{(s)}}\big(\mathbf{R^{(s)}} + \mathbf{\Sigma}^{(4)}\big)^{-1}\mathbf{y}
\qquad\text{(Sionna: }\hat{\mathbf{h}}^{(5)} = \mathbf{C}\hat{\mathbf{h}}^{(4)}\text{, }\mathbf{C}=\mathbf{R^{(s)}}\big(\mathbf{R^{(s)}}+\mathbf{\Sigma}^{(4)}\big)^{-1}\text{)}
$$

When **Σ**$^{(4)}$ is diagonal, $\mathbf{C}$ shrinks each antenna component toward the spatially correlated mean dictated by **R⁽ˢ⁾**.

#### Error variance after spatial LMMSE

$$
\mathbf{M} = \mathbf{R^{(s)}} - \mathbf{C}\,\mathbf{R^{(s)}}, \qquad
\mathbf{\Sigma}^{(5)} = \mathrm{diag}(\mathbf{M})
$$

No further scaling after the last pass.

#### Intuition

Antennas see the same scatterers → spatial correlation in **R⁽ˢ⁾**. LMMSE averaging across antennas reduces noise while preserving structure; a per-antenna LS estimate ignores cross-antenna correlation.

---

### Domain summary (`order="f-t-s"`)

| Domain | Pass | **x** (estimate) | **y** (observation) | $\mathbf{R}_{\mathbf{x}\mathbf{x}}$ | Sparse **G**? |
|--------|------|------------------|---------------------|-------------------------------------|---------------|
| Frequency | 1 | $M$ subcarriers / row | $K_n$ pilot LS / row | **R⁽ᶠ⁾** | yes ($\mathbf{\Pi}_n^{\intercal}$) |
| — | 2 | scale (freq → time) | — | **R⁽ᶠ⁾** | — |
| Time | 3 | $N$ symbols / column | $L_m$ DMRS values / column | **R⁽ᵗ⁾** | yes ($\tilde{\mathbf{\Pi}}_m^{\intercal}$) |
| — | 4 | scale (time → space) | — | **R⁽ᵗ⁾** | — |
| Spatial | 5 | $L$ antennas / RE | $L$ pass-4 estimates | **R⁽ˢ⁾** | no ($\mathbf{I}_L$) |

Our script uses `order="f-t"` only (frequency + time; passes 1–3).

---

## References

- [Sionna `LMMSEInterpolator` API](https://nvlabs.github.io/sionna/phy/api/ofdm/sionna.phy.ofdm.LMMSEInterpolator.html)
- [Sionna `tdl_freq_cov_mat`](https://nvlabs.github.io/sionna/phy/api/ofdm/sionna.phy.ofdm.tdl_freq_cov_mat.html)
- [Sionna `tdl_time_cov_mat`](https://nvlabs.github.io/sionna/phy/api/ofdm/sionna.phy.ofdm.tdl_time_cov_mat.html)
- Example: `linear_vs_lmmse_interpolator.py`, results in `readme.md`
