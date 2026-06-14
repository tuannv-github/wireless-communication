# Beamforming Demo (MATLAB)

Ray-tracing beamforming simulation on the **SUTD campus** map. A base station with an 8×8 NR panel array tracks a moving UE along a path; at each stop, rays are traced, a custom CDL channel is built, and **SVD-based beamforming weights** steer the BS and UE antenna patterns toward the strongest path.

**Demo video:** [Beamforming Demo using Matlab](https://www.youtube.com/watch?v=dybRBoWt04w)

## What it does

1. Opens a 3D **site viewer** with OpenStreetMap and building geometry from `SUTD.osm`.
2. Places a **6 GHz** base station (8×8 cross-polarized panel) and a **2×2** UE array on the campus.
3. Moves the UE along **20 interpolated positions** between a start and stop point.
4. For each position, runs **SBR ray tracing** (up to 10 reflections), builds a custom `nrCDLChannel`, and estimates the channel.
5. Computes transmit/receive beamforming weights with **`getBeamformingWeights`** (SVD on the channel matrix).
6. Updates **radiation patterns** and **ray paths** in the viewer as the UE moves.

| Parameter | Value |
|-----------|-------|
| Carrier frequency | 6 GHz |
| Bandwidth | 10 MHz (52 RB, 15 kHz SCS) |
| BS array | 8×8 NR panel, azimuth −150°, elevation −10° |
| UE array | 2×2 isotropic panel |
| Ray tracing | SBR, max 10 reflections |

## Requirements

- MATLAB R2021b or later (recommended)
- **Antenna Toolbox** — `siteviewer`, `txsite`, `rxsite`, `raytrace`, phased arrays
- **5G Toolbox** or **Communications Toolbox** — `nrCDLChannel`, `nrPerfectChannelEstimate`, NR OFDM info

## Run

Open and run `BeamForming.mlx` from this folder (keep `SUTD.osm` in the same directory):

```matlab
BeamForming
```

Run sections top to bottom. The live script animates ray paths and antenna patterns in the site viewer.

## Files

| File | Description |
|------|-------------|
| `BeamForming.mlx` | Main live script |
| `SUTD.osm` | OpenStreetMap building data for SUTD campus |
| `azemuth_0_elevation_0_Screenshot *.png` | Example viewer / pattern screenshots |

## Screenshots

![Site viewer and ray tracing](azemuth_0_elevation_0_Screenshot%202025-03-12%20135353.png)

![Antenna radiation patterns](azemuth_0_elevation_0_Screenshot%202025-03-12%20135406.png)
