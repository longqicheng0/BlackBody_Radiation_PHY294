# Step 2: Complete Uncertainty Propagation Implementation

## Overview
Full uncertainty analysis has been integrated into the Step 2 Wien's displacement law analysis. All measured quantities now include rigorously propagated uncertainties, displayed on plots with error bars, and saved to CSV with dedicated uncertainty columns.

---

## 1. Uncertainty Sources and Values

### Angle Measurement Uncertainties
- **Instrument resolution**: `U_THETA_INSTRUMENT = 0.2°` (prism spectrometer angular resolution)
- **Peak-picking uncertainty**: `U_THETA_PEAK_PICKING = 0.1052°` (from Step 1 ensemble analysis)
- **Total**: $u(\theta) = \sqrt{0.2^2 + 0.1052^2} = 0.226°$

### Voltage and Current
- **Voltage**: `U_VOLTAGE = 0.05 V` (5% typical for DMM)
- **Current**: `U_CURRENT = 0.01 A` (1% typical for DMM)

### Reference Constants
- **Room temperature**: `U_T0 = 2.0 K` (uncertainty in reference temp)
- **Reference resistance**: `U_R0 = 0.05 Ω` (5% tolerance)
- **Temperature coefficient**: `U_ALPHA0 = 0.2e-3 K⁻¹` (coeff uncertainty)

---

## 2. Uncertainty Propagation Functions

### `propagate_theta_peak_uncertainty()`
Combines instrument and peak-picking uncertainties in quadrature:
$$u(\theta) = \sqrt{u_{\text{instrument}}^2 + u_{\text{peak-picking}}^2}$$

### `propagate_lambda_uncertainty(delta_theta, u_theta)`
Propagates angle uncertainty to wavelength using finite-difference derivative:
$$u(\lambda) = \left|\frac{d\lambda}{d\theta}\right| \cdot u(\theta)$$

The derivative is computed numerically from Wien's prism equation.

### `propagate_T_uncertainty(V, I, ...)`
Propagates V, I, and reference constant uncertainties to temperature using partial derivatives.

Temperature model: $T = T_0 + \frac{(V/I)/R_0 - 1}{\alpha_0}$

Partial derivatives (computed analytically):
$$\frac{\partial T}{\partial V} = \frac{1}{R_0 \alpha_0 I}$$
$$\frac{\partial T}{\partial I} = -\frac{V}{R_0 \alpha_0 I^2}$$
$$\frac{\partial T}{\partial T_0} = 1$$
$$\frac{\partial T}{\partial R_0} = -\frac{V/I}{R_0^2 \alpha_0}$$
$$\frac{\partial T}{\partial \alpha_0} = -\frac{V/I/R_0 - 1}{\alpha_0^2}$$

Combined in quadrature:
$$u(T) = \sqrt{\sum_i \left(\frac{\partial T}{\partial x_i} \cdot u(x_i)\right)^2}$$

### `propagate_inv_T_uncertainty(T, u_T)`
Inverse temperature uncertainty:
$$u(1/T) = \frac{u(T)}{T^2}$$

### `propagate_b_uncertainty(lambda, u_lambda, T, u_T)`
Wien constant uncertainty from wavelength and temperature:
$$u(b) = \sqrt{(T \cdot u_\lambda)^2 + (\lambda \cdot u_T)^2}$$

where $\lambda$ is in meters and $b = \lambda \times T$ in m·K.

---

## 3. CSV Output Structure

New columns added to `step2_wien_results.csv`:

| Column | Unit | Description |
|--------|------|-------------|
| `voltage_V` | V | Applied voltage |
| `current_A` | A | Filament current |
| `theta_peak_deg` | ° | Detected spectrum peak angle |
| `u_theta_peak_deg` | ° | **Angle uncertainty** |
| `delta_theta_deg` | ° | Angular separation ($\theta_{\text{init}} - \theta_{\text{peak}}$) |
| `lambda_max_nm` | nm | Peak wavelength |
| `u_lambda_max_nm` | nm | **Wavelength uncertainty** |
| `T_K` | K | Computed filament temperature |
| `u_T_K` | K | **Temperature uncertainty** |
| `inv_T` | K⁻¹ | Inverse temperature |
| `u_inv_T` | K⁻¹ | **Inverse temperature uncertainty** |
| `b_mK` | m·K | Wien constant |
| `u_b_mK` | m·K | **Wien constant uncertainty** |
| `notes` | — | Analysis status |

**Example row:**
```
...,T_K=3346.78,u_T_K=209.02,inv_T=0.000299,u_inv_T=1.87e-05,b_mK=0.003000,u_b_mK=0.000263,...
```

---

## 4. Updated Plots with Error Bars

### Plot 1: `step2_lambda_vs_inv_t.png` (2-panel layout)

**Top panel: λ_max vs 1/T with weighted regression**
- **Data points**: Blue dots with ±1σ error bars (both axes)
  - Horizontal error bars: $u(1/T)$
  - Vertical error bars: $u(\lambda)$
- **Weighted fit line**: Red dashed line
  - Weights: $w_i = 1/u_{\lambda,i}^2$ (inverse variance weighting)
  - Fitted parameters with uncertainties in legend:
    - Slope: $(m \pm u_m)$ nm·K
    - Intercept: $(c \pm u_c)$ nm
- **Chi-squared statistic** (red box, top-right):
  - $\chi^2 = \sum w_i (y_i - \hat{y}_i)^2$
  - Displayed as: `χ² = XX.XX`

**Bottom panel: Residuals with error bars**
- Residuals: $y_i - \hat{y}_i$ (data minus fit)
- Error bars: $u(\lambda_i)$ (data uncertainties)
- Reference line: $y = 0$ (red dashed)
- Shaded band: ±1σ of residuals (orange region)

### Plot 2: `step2_wien_constant_vs_t.png`

- **Data points**: Green dots with ±1σ error bars (vertical)
  - Error bars represent $u(b_i)$ for each measurement
- **Literature reference**: Red dashed line at $b_0 = 2.898 \times 10^{-3}$ m·K
- **Measurement mean**: Orange dash-dot line at $\bar{b}$
- **Uncertainty band**: Orange shaded region ±1σ around mean
- **Axis formatting**: Scientific notation for y-axis (m·K values)
- **Title**: "Wien Constant vs Temperature (Testing Constancy)"

---

## 5. Console Output

Each measurement now displays uncertainties:
```
Processing: step2_10v.txt
  → V=10.0V, I=0.617A, T=3347K, λ=896.3nm, b=3.000×10⁻³ m·K
```

Summary statistics include uncertainty information:
```
======================================================================
STEP 2 ANALYSIS SUMMARY: Wien's Displacement Law
======================================================================
Files processed: 9
Successful analyses: 9

Wien constant mean:              2.4927×10⁻³ m·K
Wien constant std dev:           0.3537×10⁻³ m·K
Literature b₀:                   2.8980×10⁻³ m·K
Percent difference:              -13.99%
======================================================================
```

---

## 6. Key Features

✅ **Full error propagation** through all analysis steps  
✅ **Weighted regression** using measurement uncertainties  
✅ **Error bars on all plots** for both axes where applicable  
✅ **Fit uncertainties** (slope ± uncertainty, intercept ± uncertainty)  
✅ **Chi-squared statistics** for model validation  
✅ **Residuals analysis** with uncertainty bands  
✅ **CSV export** with all uncertainty columns  
✅ **High-resolution output** (300 dpi PNG)  
✅ **Proper axis labels** with units: λ_max (m), 1/T (K⁻¹), b (m·K)

---

## 7. Usage

Run the analysis as before:
```bash
python3 scripts/step2_wien.py
```

All uncertainties are automatically computed and included in:
- `outputs/step2_wien_results.csv` (13 data columns + 4 uncertainty columns)
- `outputs/step2_lambda_vs_inv_t.png` (2-panel plot with residuals)
- `outputs/step2_wien_constant_vs_t.png` (scatter with error bars)

---

## 8. Example: Uncertainty for One Measurement

For `step2_10v.txt` (V=10.0V, I=0.617A):

| Quantity | Value | Uncertainty | Relative |
|----------|-------|-------------|----------|
| θ_peak | 21.12° | ±0.226° | 1.07% |
| λ_max | 896.3 nm | ±55.1 nm | 6.1% |
| T | 3347 K | ±209 K | 6.2% |
| 1/T | 0.000299 K⁻¹ | ±1.87×10⁻⁵ K⁻¹ | 6.2% |
| b | 0.003000 m·K | ±0.000263 m·K | 8.8% |

The Wien constant uncertainty combines wavelength (6.1%) and temperature (6.2%) uncertainties → 8.8%.

---

## Notes

- **Weighted regression**: Uses inverse variance weighting for proper statistical analysis
- **Chi-squared**: Enables assessment of whether uncertainties are realistic (χ²/dof ≈ 1 is good)
- **Residuals plot**: Shows systematic errors and validates the linear model
- **Finite difference derivative**: Numerically stable for complex prism equations
- All plots use **matplotlib only** (no seaborn or other dependencies)
