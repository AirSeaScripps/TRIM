# Flooding

Tools for estimating coastal flood hazard using a copula-based joint probability model. The primary output is a 100-year Total Water Level (TWL) return level that correctly accounts for the statistical dependence between tides, storm surge, and wave runup.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Prerequisites](#prerequisites)
- [Setup](#setup)
- [Workflow](#workflow)
- [Script Reference](#script-reference)
- [Method](#method)
- [Configuration Parameters](#configuration-parameters)
- [Outputs](#outputs)

---

## Overview

Total Water Level at the shoreline has three primary components:

$$\text{TWL} = \eta_{\text{tide}} + \eta_{\text{surge}} + R_2$$

where $\eta_{\text{tide}}$ is the astronomic tide, $\eta_{\text{surge}}$ is storm surge (the non-tidal residual after removing sea level rise), and $R_2$ is the 2% exceedance wave runup level from [Stockdon et al. (2006)](#stockdon-2006).

A naive approach to estimating the 100-year TWL combines separately-derived return levels for each component, which implicitly assumes they are statistically independent. In practice, large storms simultaneously produce elevated surge **and** large waves, so this assumption can significantly undercount joint occurrence probability. This toolbox instead fits a **multivariate copula** to concurrent observations of all three variables and draws from their joint distribution in a Monte Carlo simulation, giving a physically more correct extreme value estimate.

---

## Repository Structure

```
Flooding/
├── TWL_Analysis_copula.m       # Main analysis script (start here)
├── tide_data.m                 # Downloads NOAA tide/water-level data
├── Get_ERA5_Waves.m            # Downloads ERA5 wave reanalysis (calls Python)
├── Get_ERA5_Waves.py           # Python CDS API helper (called by Get_ERA5_Waves.m)
├── Stockdon2006.m              # Wave runup parameterisation (sandy beach)
├── Becker2026.m                # Wave runup parameterisation (dissipative fringing reef)
├── ddm2dd.m                    # Degrees-decimal-minutes to decimal-degrees coordinate parser
├── navd88_to_egm2008.m         # NAVD88 -> EGM2008 datum conversion via NOAA VDatum API
├── Analyze_FloodOutput.m       # Post-processing / coastline elevation analysis
├── Data/
│   └── CACC_CoastlineElevations.mat
└── WaveFlood/                  # Python virtual environment (not tracked by git)
```

---

## Prerequisites

### MATLAB

- **MATLAB R2020b or later** (earlier releases may work but are untested)
- **Statistics and Machine Learning Toolbox** — required for `copulafit`, `copularnd`, `fitdist`, `icdf`, `tiedrank`, `poissrnd`, `tcdf`, `prctile`
- **Signal Processing Toolbox** — required for the Butterworth low-pass filter (`butter`, `zp2sos`, `filtfilt`)

### Python

Wave data retrieval uses the [Copernicus Climate Data Store (CDS) API](https://cds.climate.copernicus.eu/) via a Python helper script. You need:

- **Python 3.9** (the virtual environment `WaveFlood/` is built against 3.9; other versions may require a fresh venv — see [Setup](#setup))
- **`cdsapi` package** (installed into `WaveFlood/`)
- A valid **CDS API key** stored in `~/.cdsapirc`

If NOAA tide data is sufficient and you are not computing wave runup (`CALC_WAVES = 0`), the Python environment is not needed.

---

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/AirSeaScripps/Flooding.git
cd Flooding
```

### 2. Create the Python virtual environment

```bash
python3.9 -m venv WaveFlood
source WaveFlood/bin/activate      # macOS / Linux
# WaveFlood\Scripts\activate       # Windows
pip install cdsapi
```

> **Note:** The `WaveFlood/` directory is intentionally excluded from version control (`.gitignore`). You must create it locally before running wave analysis.

### 3. Configure the CDS API key

Register at [https://cds.climate.copernicus.eu](https://cds.climate.copernicus.eu) and create `~/.cdsapirc`:

```
url: https://cds.climate.copernicus.eu/api
key: <your-uid>:<your-api-key>
```

Full instructions are on the [CDS API how-to page](https://cds.climate.copernicus.eu/how-to-api).

### 4. Python interpreter path

`Get_ERA5_Waves.m` points MATLAB at the `WaveFlood/` virtual environment automatically, resolved relative to the script's own location — no path editing required as long as `WaveFlood/` lives at the repo root as created in step 2.

---

## Workflow

```
NOAA Tides & Currents API          ERA5 / CDS API
        │                                  │
        ▼                                  ▼
   tide_data.m                   Get_ERA5_Waves.m
   (observed + predicted tides)  (Hs, Tp, wave direction)
        │                                  │
        └──────────────┬───────────────────┘
                       ▼
           TWL_Analysis_copula.m
           ┌──────────────────────────────┐
           │  1. Low-pass filter → SLR    │
           │  2. POT / GPD wave analysis  │
           │  3. Fit copula               │
           │  4. Monte Carlo (100 yr)     │
           │  5. GEV fit → 100-yr TWL     │
           └──────────────────────────────┘
```

**Quickstart:**

1. Open `TWL_Analysis_copula.m` in MATLAB.
2. Fill in the **USER INPUT** section (station ID, datum values, lat/lon, beach slope).
3. Set `CALC_WAVES = 1` if offshore wave data are available, `0` otherwise.
4. Run the script. Data are downloaded automatically from NOAA and ERA5; figures and printed results are generated at the end.

---

## Script Reference

### `TWL_Analysis_copula.m`

The main analysis script. Orchestrates data download, signal decomposition, extreme value analysis, copula fitting, and Monte Carlo simulation.

**Inputs** (set in the USER INPUT block):

| Variable | Description |
|---|---|
| `USE_NOAA` | `1` = download NOAA tide gauge data; `0` = skip the NOAA download |
| `years` | Year range for NOAA download, e.g. `[1955:2025]'` |
| `station` | NOAA station ID string, e.g. `'9410230'` |
| `MHHW` | Mean Higher High Water height above NAVD88 (m) |
| `NAVD88` | Additional NAVD88 offset if datum is not directly available (m) |
| `NAVD88_to_EGM` | Offset from NAVD88 to EGM2008 geoid at the site (m) — see `navd88_to_egm2008.m` |
| `CALC_WAVES` | `1` = trivariate copula with wave runup; `0` = bivariate, no waves |
| `WAVE_MODEL` | `'Stockdon'` (sandy beach, default) or `'Becker'` (dissipative fringing reef) — see `Stockdon2006.m` / `Becker2026.m` |
| `beta` | Beach slope for the Stockdon runup model (used when `WAVE_MODEL = 'Stockdon'`) |
| `D0` | Mean reef-flat bed elevation relative to the tide gauge datum (used when `WAVE_MODEL = 'Becker'`) |
| `lat`, `lon` | ERA5 wave extraction point |
| `Hs_thresh_pct` | Percentile threshold for POT wave analysis (default: 95) |
| `Hs_max` | Physical upper bound on sampled $H_s$ in the Monte Carlo (m); set to `Inf` to disable |
| `COPULA_FAMILY` | `'t'` (Student-t, recommended) or `'Gumbel'` (bivariate only) |
| `USE_PRODUCT_COPULA` | `false` (default) fits the copula normally; `true` replaces it with independent draws, for comparison against the independence assumption |
| `repeats` | Number of 100-year Monte Carlo repetitions (default: 300) |
| `tau_hours` | Inter-storm separation window for POT declustering (hours) |

---

### `tide_data.m`

Downloads hourly observed and predicted water levels from the [NOAA CO-OPS API](https://tidesandcurrents.noaa.gov/api-helper/url-generator.html), referenced to MHHW datum.

```matlab
[predi, obsv, d_t] = tide_data('0101', '1231', years, station)
```

| Argument | Description |
|---|---|
| `begin_date` | Start of each annual window as `'MMDD'` |
| `end_date` | End of each annual window as `'MMDD'` |
| `years` | Column vector of years |
| `station` | NOAA station ID string |

**Returns:** `predi` (predicted tide, m), `obsv` (observed water level, m), `d_t` (datetime vector).

---

### `Get_ERA5_Waves.m` + `Get_ERA5_Waves.py`

Downloads the most recent 30 full calendar years of ERA5 offshore wave data from the Copernicus CDS, computed dynamically from the current date (e.g. run in 2038 → 2008–2037). The MATLAB function sets up the Python environment and iterates over three variables: significant wave height ($H_s$), peak wave period ($T_p$), and mean wave direction.

```matlab
[Hs_ERA5, Tp_ERA5, time_ERA5, D_ERA5] = Get_ERA5_Waves(lat, lon)
```

The function requests a 1° bounding box around the target point to guard against the target location falling on a land cell, then selects the nearest valid ocean grid point.

> ERA5 data are available at 00:00 and 12:00 UTC, giving a 12-hourly time series well suited for copula fitting at storm timescales.

---

### `Stockdon2006.m`

Vectorised implementation of the [Stockdon et al. (2006)](https://doi.org/10.1016/j.coastaleng.2005.12.005) empirical runup model.

```matlab
[R2, setup, swash, sinc, sig] = Stockdon2006(Hs, Tp, beta)
```

For reflective/intermediate beaches ($\zeta \geq 0.3$):

$$R_2 = 1.1 \left( 0.35\,\beta\sqrt{H_s L_0} + \frac{\sqrt{H_s L_0 (0.563\,\beta^2 + 0.004)}}{2} \right)$$

For dissipative conditions ($\zeta < 0.3$):

$$R_2 = 0.043\sqrt{H_s L_0}$$

where the deep-water wavelength is $L_0 = \dfrac{g}{2\pi}T_p^2$ and the Iribarren number is $\zeta = \beta / \sqrt{H_s / L_0}$.

---

### `Becker2026.m`

Alternative wave-driven shoreline water level model for dissipative fringing reefs, following Becker & Merrifield (2026). Drop-in replacement for `Stockdon2006.m` when `WAVE_MODEL = 'Becker'` — validated at Ipan reef, Guam, and not appropriate for wave-like reefs or non-reef (sandy beach / open coast) sites.

```matlab
[eta2, eta_bar, eta_prime, hr] = Becker2026(H0, Tp, tide, surge, D0)
```

---

### `ddm2dd.m`

Parses degrees-decimal-minutes coordinate strings (as published on NOAA tide gauge pages, e.g. `"48° 6.7 N"`) into signed decimal degrees for use as `lat`/`lon` inputs elsewhere in the pipeline.

---

### `navd88_to_egm2008.m`

Converts NAVD88 orthometric heights to EGM2008 (WGS84) ellipsoidal heights via NOAA's VDatum REST API — used to compute `NAVD88_to_EGM` without manually running points through the VDatum web GUI.

```matlab
[zOut, info] = navd88_to_egm2008(lon, lat, zNAVD88, 'region', 'westcoast')
```

---

### `Analyze_FloodOutput.m`

Post-processing script for flood output analysis using CACC coastline elevation data (`Data/CACC_CoastlineElevations.mat`). Work in progress.

---

## Method

### Signal decomposition

The observed water level $\eta_{\text{obs}}$ is decomposed as:

$$\eta_{\text{obs}} = \eta_{\text{tide}} + \eta_{\text{SL}} + \eta_{\text{surge}}$$

where $\eta_{\text{tide}}$ is the NOAA-predicted tidal signal, $\eta_{\text{SL}}$ is a low-frequency sea level rise trend (extracted with a 4th-order Butterworth low-pass filter with a 12-month cutoff), and $\eta_{\text{surge}} = \eta_{\text{obs}} - \eta_{\text{tide}} - \eta_{\text{SL}}$ is the storm surge residual.

### Wave extreme value analysis (POT/GPD)

Wave extremes are estimated using Peaks Over Threshold (POT). Consecutive $H_s$ exceedances above the $p$-th percentile threshold $u$ are declustered with a minimum inter-storm gap of `tau_hours`, keeping only the peak $H_s$ of each storm cluster. Storm-peak excesses $e = H_s - u$ are fitted with a Generalised Pareto Distribution (GPD):

$$F(e;\,\sigma,\xi) = 1 - \left(1 + \frac{\xi\, e}{\sigma}\right)^{-1/\xi}$$

The Poisson-GPD 100-year return level is:

$$H_{s,100} = u + \text{GPD}^{-1}\!\left(1 - \frac{1}{\lambda\, T}\right)$$

where $\lambda$ is the storm rate (storms yr⁻¹) and $T = 100$ yr. **This value is a diagnostic only** — the 100-year TWL is computed directly inside the Monte Carlo below.

### Copula fitting

The code supports two modes:

**Trivariate (preferred, `CALC_WAVES = 1`):** A Student-t copula is fitted to concurrent 12-hourly triplicates $(\eta_{\text{tide}},\, \eta_{\text{surge}},\, H_s)$ after transforming each variable to uniform marginals via rank-based pseudo-observations (Weibull plotting position). The t-copula is parameterised by a $3 \times 3$ correlation matrix $\boldsymbol{\rho}$ and degrees of freedom $\nu$:

$$C(\mathbf{u};\,\boldsymbol{\rho},\nu) = t_{\boldsymbol{\rho},\nu}\!\left(t_\nu^{-1}(u_1),\, t_\nu^{-1}(u_2),\, t_\nu^{-1}(u_3)\right)$$

As $\nu \to \infty$ the t-copula recovers the Gaussian (no tail dependence) limit. Finite $\nu$ produces symmetric upper and lower tail dependence with coefficient:

$$\lambda_U = 2\,t_{\nu+1}\!\left(-\sqrt{\frac{(\nu+1)(1-\rho)}{1+\rho}}\right)$$

**Bivariate (fallback, `CALC_WAVES = 0`):** A bivariate t- or Gumbel copula is fitted to $(\eta_{\text{tide}},\, \eta_{\text{surge}})$ only. Wave runup is set to zero.

### $T_p$ regression

Rather than fixing $T_p$ at a single design value, a power-law regression is fitted to storm-peak $(H_s, T_p)$ pairs in log-space:

$$\ln T_p = a\,\ln H_s + b \quad \Longrightarrow \quad T_p = e^b\,H_s^a$$

Log-space residuals $\varepsilon \sim \mathcal{N}(0,\sigma_{T_p}^2)$ are sampled independently for each Monte Carlo draw to propagate $T_p$ uncertainty into the runup calculation.

### Monte Carlo simulation

For each of `repeats` repetitions:

1. Draw storm count $N_y \sim \text{Poisson}(\lambda)$ for each of 100 simulated years.
2. Sample $N_{\text{total}} = \sum_y N_y$ correlated uniform triplicates from the fitted copula.
3. Back-transform to physical units using empirical quantile functions for tide and surge, and a hybrid empirical/GPD quantile function for $H_s$.
4. Compute $T_p$ from the regression (with lognormal scatter) and $R_2$ from Stockdon (2006).
5. Sum to get $\text{TWL} = \eta_{\text{tide}} + \eta_{\text{surge}} + R_2$; record the annual maximum.

Annual maxima across all repetitions are pooled and fitted with a Generalised Extreme Value (GEV) distribution:

$$F(x;\,\mu,\sigma,\xi) = \exp\!\left[-\left(1 + \xi\,\frac{x - \mu}{\sigma}\right)^{-1/\xi}\right]$$

The **100-year TWL return level** is $\text{GEV}^{-1}(0.99)$, reported relative to MHHW and EGM2008.

---

## Configuration Parameters

The parameters most likely to need adjustment for a new site are:

| Parameter | Where | Notes |
|---|---|---|
| `station` | USER INPUT | NOAA station ID — find yours at [tidesandcurrents.noaa.gov](https://tidesandcurrents.noaa.gov/map/index.html) |
| `years` | USER INPUT | Use at least 30 years for reliable statistics |
| `MHHW` | USER INPUT | MHHW height above NAVD88 in metres — available from the NOAA station datums page |
| `NAVD88_to_EGM` | USER INPUT | Compute at [vdatum.noaa.gov](https://vdatum.noaa.gov/vdatumweb/) |
| `lat`, `lon` | USER INPUT | Offshore point for ERA5 wave extraction |
| `beta` | USER INPUT | Beach slope — estimate from lidar/DEM; typical sandy beach range 0.03–0.15 |
| `Hs_thresh_pct` | USER INPUT | POT threshold percentile; 90–98 is a reasonable range (default 95) |
| `repeats` | USER INPUT | Higher values reduce Monte Carlo variance; 300 is a practical default |

---

## Outputs

The script prints a results summary to the MATLAB console, including:

- **100-year TWL** relative to MHHW and EGM2008 (m)
- Fitted copula parameters ($\boldsymbol{\rho}$, $\nu$)
- Pairwise upper tail dependence coefficients $\lambda_U$ with a plain-language classification (near-independent / weak / IMPORTANT)
- Sea level rise rate (mm yr⁻¹)
- Storm rate $\lambda$ and diagnostic $H_{s,100}$, $T_{p,100}$

Figures generated (all rendered in the MATLAB figure window):

1. GPD fit to storm-peak $H_s$ excesses
2. $H_s$–$T_p$ scatter with fitted power-law regression and uncertainty envelopes
3. Pairwise copula scatter plots (uniform marginals)
4. Chi-plots of upper tail dependence for each variable pair
5. GEV fit to simulated TWL annual maxima
6. Sea level rise time series with linear trend
7. Storm surge scatter
8. Annual maxima of observed water level (NOAA-method comparison)

---

## References

**Stockdon 2006**
Stockdon, H. F., Holman, R. A., Howd, P. A., & Sallenger, A. H. (2006). Empirical parameterization of setup, swash, and runup. *Coastal Engineering*, 53(7), 573–588. https://doi.org/10.1016/j.coastaleng.2005.12.005

**ERA5**
Hersbach, H., et al. (2020). The ERA5 global reanalysis. *Quarterly Journal of the Royal Meteorological Society*, 146, 1999–2049. https://doi.org/10.1002/qj.3803
