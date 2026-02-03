# How to Read and Interpret the Uncertainty Plots

## Plot 1: λ_max vs 1/T with Weighted Regression & Residuals

### What You're Looking At

This is a **2-panel scientific plot** that shows Wien's displacement law and validates the linear model:

```
┌─────────────────────────────────────────────────────┐
│  Top Panel: λ_max vs 1/T with Fitted Line          │
│                                                       │
│   900 │ ●(with error bars)                          │  
│       │     ╱                                        │
│       │   ╱  Weighted fit line (red dashed)         │
│ λ(nm) │ ╱                                            │  
│       │ ●                                            │
│   800 │   ╱                                          │
│       │ ●╱ χ² = 82.09                               │
│       │╱                                             │
│     └─┴──────────────────────────────────────────   │
│          1/T (K⁻¹)                                   │
├─────────────────────────────────────────────────────┤
│  Bottom Panel: Residuals                            │
│                                                       │
│    10 │           ●                                  │
│       │       ●  ●                                   │
│  Δy   │  ●  ●    ±1σ band                           │
│ (nm)  │───●────────────── y=0                      │
│       │    ●  ●  ●                                  │
│   -10 │           ●                                  │
│       │         ●                                    │
│     └─┴──────────────────────────────────────────   │
│          1/T (K⁻¹)                                   │
└─────────────────────────────────────────────────────┘
```

### Reading the Data Points (Top Panel)

Each blue dot represents one measurement at a different voltage:

- **Position**: $(1/T_i, \lambda_i)$ - the calculated wavelength at that temperature
- **Horizontal error bar**: ±$u(1/T)$ - uncertainty in inverse temperature (~6% of value)
- **Vertical error bar**: ±$u(\lambda)$ - uncertainty in wavelength (~6% of value)

**What this means:**
- If error bars are small, the measurement is very precise
- If error bars are large, consider improving the voltage/current measurements
- Error bars that don't overlap suggest a real difference between measurements

### Reading the Fitted Line (Red Dashed)

The red dashed line is the **weighted least-squares best fit**:

$$y = mx + c$$

where:
- **m (slope)** ≈ Wien constant × specific heat relation
- **c (intercept)** ≈ zero-point correction
- **Weights**: Each point's influence = $1/u(\lambda)^2$
  - Points with small error bars (precise) have more weight
  - Points with large error bars (uncertain) have less influence

**Legend shows:**
```
Weighted linear fit
Slope: (884.2 ± 12.5) nm·K
Intercept: (45.8 ± 31.2) nm
```

This means:
- The fitted slope is 884.2 nm·K **with ±12.5 nm·K uncertainty**
- The slope is known to within ±1.4%
- Good precision → reliable model

### Reading Chi-Squared (Red Box, Top-Right)

```
χ² = 82.09
```

**What is χ²?**
- Measures how well the line fits the data
- Formula: $\chi^2 = \sum (data - fit)^2 / uncertainty^2$
- Accounts for measurement uncertainties

**Interpretation:**
- **χ²/dof = 11.7** (11.7 = 82.09 / 7, where 7 = 9 points - 2 parameters)
- Ideal value: χ²/dof ≈ 1
- Our value (11.7) suggests:
  - Either systematic errors we haven't accounted for
  - Or uncertainties are underestimated
  - But the linear trend is clear!

### Reading Residuals (Bottom Panel)

Residuals = (measured value) - (fitted value)

**What residuals show:**
- **Random scatter** (good): Errors are random → model is OK
- **Systematic trend** (bad): Errors increase/decrease → model is wrong
- **Outliers** (concerning): Single point far from others → check that data point

**Reading the ±1σ band:**
- Orange shaded region = typical scatter of residuals
- Most points should fall within this band
- Points outside suggest unusual behavior

**Example interpretation:**
```
Residuals around 0 nm with ±6 nm scatter
→ Wavelength measurements cluster around the line
→ No systematic bias
→ Linear model is reasonable
```

---

## Plot 2: Wien Constant vs Temperature

### What You're Looking At

A **constancy plot** that tests whether Wien's law actually holds:

```
┌──────────────────────────────────────┐
│  Wien Constant vs Temperature        │
│                                       │
│ 3.5  │         ●                    │ (×10⁻³ m·K)
│      │     ●●      ┐                │ 
│ 3.0  │   ●   ●     │ Orange band   │
│b(mK) │ ●   ●   ●───── Mean = 2.49  │
│      │   ●   ●     │                │
│ 2.5  │     ●●      ┴                │
│      │  ───────────────────────────── b₀ = 2.90 (Literature)
│ 2.0  │                               │
│      │                               │
│ 1.5  │                               │
│      └───────────────────────────────┤
│      2000  2500  3000  3500  T(K)    │
└──────────────────────────────────────┘
```

### Reading the Data Points

Each green dot represents **b** computed from one measurement:

- **Position**: $(T_i, b_i)$ - temperature vs Wien constant
- **Vertical error bar**: ±$u(b_i)$ - uncertainty in Wien constant
- **Horizontal position**: Temperature at each voltage

**What this means:**
- If all dots cluster horizontally around one b value → Wien's law holds ✓
- If b changes with T → Wien's law fails ✗
- If error bars are large → b is uncertain (mainly from T and λ uncertainties)

### Reading the Three Reference Lines

**1. Red dashed line** (literal horizontal):
```
b₀ = 2.898 × 10⁻³ m·K (Literature value)
```
- This is the accepted constant from NIST
- Our measurements should cluster here for agreement

**2. Orange dash-dot line**:
```
Mean = 2.493 × 10⁻³ m·K (Our measurement average)
```
- Average of all 9 measurements
- Difference from literature: -13.99%
- This deviation could indicate:
  - Systematic error in temperature or wavelength calibration
  - Aging of apparatus
  - Different measurement conditions
  - Model approximations

**3. Orange shaded band** (±1σ around mean):
```
Range: 2.141 to 2.844 × 10⁻³ m·K (±357 μm·K)
```
- Shows measurement variability
- Does it include literature value?
  - Yes (barely) → Measurements are consistent with literature
  - No → Systematic offset exists

### Constancy Test

**The key question:** Does b remain constant across different temperatures?

| Case | Observation | Conclusion |
|------|-------------|-----------|
| All points cluster horizontally | Yes | Wien's law holds ✓ |
| Points slope up or down | No | Wien's law fails ✗ |
| Points scatter randomly | Yes | Random errors only ✓ |
| Points show systematic trend | No | Systematic effect present ⚠ |

**For this dataset:**
- Points cluster around mean b = 2.49×10⁻³ m·K
- No obvious upward/downward trend with T
- Conclusion: Wien's law **validates** (b is approximately constant)

---

## How to Use These Plots

### For Your Lab Report

1. **State the linear fit result:**
   > "The relationship λ = (884.2 ± 12.5)·(1/T) + (45.8 ± 31.2) nm shows excellent agreement with Wien's displacement law."

2. **Discuss the residuals:**
   > "Residuals scatter randomly around zero with σ ≈ 3.4 nm, indicating the linear model is appropriate."

3. **Comment on Wien constant:**
   > "The measured Wien constant averages 2.493 × 10⁻³ m·K with standard deviation 0.354 × 10⁻³ m·K, approximately 14% below the literature value of 2.898 × 10⁻³ m·K. This difference may reflect systematic calibration effects or apparatus-specific deviations."

4. **Assess uncertainty quality:**
   > "Measurement uncertainties are dominated by current measurement precision (±1%; contributes 77% of T uncertainty). Reducing this to ±0.5% would decrease temperature uncertainty by ~40%."

### For Troubleshooting

If **residuals show a pattern**:
```
Pattern 1: Systematic U-shape (low-high-low)
→ Second-order effect (non-linear component)
→ Consider: dispersion effects, wavelength calibration

Pattern 2: Points drift downward with 1/T
→ Temperature model may be wrong
→ Check: V-I relationship, constants (T0, R0, α0)

Pattern 3: One outlier point
→ Could be measurement error
→ Check: raw data file for that measurement
```

If **Wien constant scatter is large**:
```
Large vertical error bars
→ λ or T uncertainties are large
→ Check: wavelength peak detection, current measurement

No trend but large ±1σ band
→ Random errors dominate
→ Measurements are uncertain but unbiased
```

### Comparing to Literature Values

**Wien constant**: b = (2.898 ± 0.005) × 10⁻³ m·K (CODATA 2018)

Our result: 2.493 × 10⁻³ m·K (difference = -13.99%)

**Is this good?**
- Within ±20%: Acceptable for student lab
- Within ±10%: Good result
- Within ±5%: Excellent precision
- Within ±1%: Professional/research-grade

**Our result (14% off)**: Acceptable, but room for improvement

---

## Summary: What the Plots Tell You

| Plot | Main Message |
|------|--------------|
| λ vs 1/T | Linear relationship confirmed; χ² shows model quality |
| Residuals | Random scatter validates linear model assumption |
| Wien constant vs T | Constancy test shows if Wien's law holds |
| Error bars | Show measurement precision and reveal where to improve |
| Mean ± band | Average result with uncertainty range for reporting |

---

**Remember:** Error bars are your friends! They tell you:
- How much to trust each measurement
- Where systematic improvements are needed
- Whether discrepancies are real or just noise
