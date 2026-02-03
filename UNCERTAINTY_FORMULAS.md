# Uncertainty Propagation Formulas - Step 2 Wien Analysis

## 1. Angle Measurement Uncertainty

### Combining Independent Sources
$$u(\theta) = \sqrt{u_{\text{instrument}}^2 + u_{\text{peak-picking}}^2}$$

$$u(\theta) = \sqrt{(0.2°)^2 + (0.1052°)^2} = 0.226°$$

---

## 2. Wavelength Uncertainty from Angle

### Using Derivative Method
$$u(\lambda) = \left|\frac{d\lambda}{d\theta}\right| \cdot u(\theta)$$

The Wien prism equation is:
$$\lambda = \sqrt{\frac{A}{\sqrt{(c_1 \sin\theta + c_2)^2 + 0.75} - B}}$$

where $A = 13900$ nm², $B = 1.689$, $c_1 = 2/\sqrt{3}$, $c_2 = 0.5$

The derivative is computed using **finite difference**:
$$\frac{d\lambda}{d\theta} \approx \frac{\lambda(\theta + h) - \lambda(\theta - h)}{2h} \quad (h = 10^{-4}°)$$

**For the 10V measurement:** $|d\lambda/d\theta| \approx 244.4$ nm/degree

$$u(\lambda_{\text{max}}) = 244.4 \text{ nm/degree} \times 0.226° = 55.1 \text{ nm}$$

---

## 3. Temperature Uncertainty (Most Complex)

### Temperature Model
$$T = T_0 + \frac{(V/I)/R_0 - 1}{\alpha_0}$$

where:
- $T_0 = 293.0$ K (room temperature reference)
- $V$ = applied voltage (V)
- $I$ = filament current (A)
- $R_0 = 1.1$ Ω (reference resistance)
- $\alpha_0 = 4.5 \times 10^{-3}$ K⁻¹ (temperature coefficient)

### Partial Derivatives
$$\frac{\partial T}{\partial V} = \frac{1}{R_0 \alpha_0 I}$$

$$\frac{\partial T}{\partial I} = -\frac{V}{R_0 \alpha_0 I^2}$$

$$\frac{\partial T}{\partial T_0} = 1$$

$$\frac{\partial T}{\partial R_0} = -\frac{V/I}{R_0^2 \alpha_0}$$

$$\frac{\partial T}{\partial \alpha_0} = -\frac{V/I/R_0 - 1}{\alpha_0^2}$$

### Combined Uncertainty (Quadrature)
$$u(T) = \sqrt{\left(\frac{\partial T}{\partial V} u_V\right)^2 + \left(\frac{\partial T}{\partial I} u_I\right)^2 + \left(\frac{\partial T}{\partial T_0} u_{T_0}\right)^2 + \left(\frac{\partial T}{\partial R_0} u_{R_0}\right)^2 + \left(\frac{\partial T}{\partial \alpha_0} u_{\alpha_0}\right)^2}$$

**For the 10V measurement** ($V=10.0$ V, $I=0.617$ A):

| Component | Value | Contribution |
|-----------|-------|--------------|
| $\frac{\partial T}{\partial V} u_V$ | 19.3 K | 0.8% |
| $\frac{\partial T}{\partial I} u_I$ | 161.7 K | 77.7% |
| $\frac{\partial T}{\partial T_0} u_{T_0}$ | 2 K | 0.1% |
| $\frac{\partial T}{\partial R_0} u_{R_0}$ | 32.8 K | 5.0% |
| $\frac{\partial T}{\partial \alpha_0} u_{\alpha_0}$ | 75.8 K | 16.4% |

$$u(T) = \sqrt{19.3^2 + 161.7^2 + 2^2 + 32.8^2 + 75.8^2} = 209.0 \text{ K}$$

**Key finding:** Current uncertainty dominates (77.7% of total variance)

---

## 4. Inverse Temperature Uncertainty

### Using Chain Rule
$$u(1/T) = \left|\frac{d(1/T)}{dT}\right| \cdot u(T) = \frac{u(T)}{T^2}$$

**For the 10V measurement:**
$$u(1/T) = \frac{209.0 \text{ K}}{(3347 \text{ K})^2} = 1.87 \times 10^{-5} \text{ K}^{-1}$$

---

## 5. Wien Constant Uncertainty

### Wien Constant Definition
$$b = \lambda \times T$$

where $\lambda$ is in meters and $T$ in Kelvin, giving $b$ in m·K.

### Uncertainty Propagation
$$u(b) = \sqrt{\left(\frac{\partial b}{\partial \lambda} u_\lambda\right)^2 + \left(\frac{\partial b}{\partial T} u_T\right)^2}$$

$$u(b) = \sqrt{(T \cdot u_\lambda)^2 + (\lambda \cdot u_T)^2}$$

**For the 10V measurement:**
- $\lambda = 896.3 \times 10^{-9}$ m = $8.963 \times 10^{-7}$ m
- $u_\lambda = 55.1 \times 10^{-9}$ m = $5.51 \times 10^{-8}$ m
- $T = 3346.78$ K
- $u_T = 209.02$ K

$$u(b) = \sqrt{(3346.78 \times 5.51 \times 10^{-8})^2 + (8.963 \times 10^{-7} \times 209.02)^2}$$

$$u(b) = \sqrt{(1.845 \times 10^{-4})^2 + (1.873 \times 10^{-4})^2}$$

$$u(b) = \sqrt{3.403 \times 10^{-8} + 3.509 \times 10^{-8}} = 2.63 \times 10^{-4} \text{ m·K}$$

**Result:** $b = (3.000 \pm 0.263) \times 10^{-3}$ m·K

---

## 6. Weighted Linear Regression

### Least Squares Problem
For data points $(x_i, y_i)$ with uncertainties $u_i$ on the y-axis:

$$\text{minimize} \quad \chi^2 = \sum_{i=1}^{n} w_i (y_i - (mx_i + c))^2$$

where weights $w_i = 1/u_i^2$

### Solution (Normal Equations)
$$\begin{bmatrix} \sum w_i x_i^2 & \sum w_i x_i \\ \sum w_i x_i & \sum w_i \end{bmatrix} \begin{bmatrix} m \\ c \end{bmatrix} = \begin{bmatrix} \sum w_i x_i y_i \\ \sum w_i y_i \end{bmatrix}$$

### Parameter Uncertainties
From the covariance matrix $\mathbf{C} = (\mathbf{X}^T \mathbf{W} \mathbf{X})^{-1}$:

$$u(m) = \sqrt{C_{11}}$$
$$u(c) = \sqrt{C_{22}}$$

### Chi-squared Goodness-of-Fit
$$\chi^2 = \sum_{i=1}^{n} w_i (y_i - \hat{y}_i)^2$$

$$\chi^2_{\text{reduced}} = \frac{\chi^2}{n - k}$$

where $k$ = number of fitted parameters (2 for linear fit), $n$ = number of data points

**For our Wien analysis** ($n=9$, $k=2$):
- $\chi^2 = 82.09$
- $\chi^2_{\text{reduced}} = 82.09 / 7 = 11.7$

$\chi^2_{\text{red}} \approx 1$ indicates good fit; $\chi^2_{\text{red}} > 3$ suggests systematic problems or underestimated uncertainties.

---

## 7. Residual Standard Deviation

### Estimate of Fit Quality
$$\sigma_{\text{residuals}} = \sqrt{\frac{\chi^2}{n-k}}$$

This represents the typical scatter of points around the best-fit line.

$$\sigma_{\text{residuals}} = \sqrt{\frac{82.09}{7}} = 3.42 \text{ nm}$$

---

## Summary: All Uncertainties for One Measurement

| Quantity | Value | Uncertainty | Method |
|----------|-------|-------------|--------|
| $\theta_{\text{peak}}$ | 21.12° | ±0.226° | Quadrature: instrument + peak-picking |
| $\lambda_{\text{max}}$ | 896.3 nm | ±55.1 nm | Derivative: $d\lambda/d\theta \times u(\theta)$ |
| $T$ | 3347 K | ±209 K | Partial derivatives + quadrature |
| $1/T$ | 0.000299 K⁻¹ | ±1.87×10⁻⁵ K⁻¹ | Chain rule: $u(1/T) = u(T)/T^2$ |
| $b$ | 0.003000 m·K | ±0.000263 m·K | Quadrature: $(Tu_\lambda)^2 + (\lambda u_T)^2$ |

---

## Comparison to Measurements

**Expected uncertainty ranges** based on instrument specs:

| Instrument | Spec | Used | Rationale |
|-----------|------|------|-----------|
| Digital Multimeter (V) | ±(0.5% + 1 digit) | ±0.05 V | Conservative estimate |
| Digital Multimeter (I) | ±(1% + 1 digit) | ±0.01 A | Conservative estimate |
| Prism spectrometer | ±0.1° typical | ±0.2° | Over-specified for safety |
| Temperature ref | ~±2 K | ±2 K | Lab ambient uncertainty |
| Resistor (1.1 Ω, 1%) | ±0.011 Ω | ±0.05 Ω | Conservative for aging |

---

**Note**: All formulas implemented in Python using numpy for numerical stability and scipy for optimization where applicable.
