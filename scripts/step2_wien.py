#!/usr/bin/env python3
"""
Blackbody Lab - Step 2 Analysis: Wien's Displacement Law

Processes Step 2 voltage-dependent scans to extract peak wavelength and compute
Wien constant, filament temperature, and validate Wien's displacement law.

Author: PHY294 Lab Analysis Tool
"""

import argparse
import re
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use('Agg')  # Non-GUI backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Optional scipy for smoothing
try:
    from scipy.signal import savgol_filter, find_peaks
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


# ============================================================================
# CONFIGURATION
# ============================================================================

class Config:
    """Configuration for Step 2 Wien analysis."""
    DATA_DIR = "Blackbody_Lab_Data"
    STEP2_DIR = "Blackbody_Lab_Data"  # Step 2 files are directly in data dir
    PATTERN = "**/step2*.txt"
    OUTPUT_DIR = "outputs"
    
    # Calibration from Step 1
    THETA_INIT_DEG = 78.28  # Unrefracted light reference angle
    
    # Peak detection
    PEAK_ANGLE_WINDOW = (10.0, 35.0)  # Search window for main peak (degrees)
    SMOOTHING_WINDOW = 11
    SMOOTHING_POLY = 3
    PROMINENCE_THRESHOLD = 0.01  # Volts
    
    # Wien's displacement law constants
    WIEN_A = 13900  # nm^2
    WIEN_B = 1.689
    
    # Temperature formula constants
    T0 = 293.0  # K
    R0 = 1.1    # ohm
    ALPHA0 = 4.5e-3  # 1/K
    
    # Wien constant literature value
    WIEN_B0 = 2.898e-3  # m·K
    
    # ========== UNCERTAINTIES ==========
    # Instrument/angle resolution
    U_THETA_INSTRUMENT = 0.2  # degrees (prism spectrometer angular resolution)
    
    # Peak-picking uncertainty (std dev from Step 1 analysis)
    U_THETA_PEAK_PICKING = 0.1052  # degrees (from Step 1 ensemble uncertainty)
    
    # Voltage and current measurement uncertainties
    U_VOLTAGE = 0.05  # V (5% or 50mV, whichever is larger - typical for DMM)
    U_CURRENT = 0.01  # A (1% or 10mA - typical for DMM)
    
    # Reference constants uncertainties
    U_T0 = 2.0  # K (room temperature uncertainty)
    U_R0 = 0.05  # ohm (5% tolerance on reference resistor)
    U_ALPHA0 = 0.2e-3  # 1/K (temperature coefficient uncertainty)
    
    DEBUG = False


# ============================================================================
# DATA LOADING
# ============================================================================

def load_scan(file_path: Path) -> pd.DataFrame:
    """
    Load a Step 2 scan file (tab-separated).
    
    Args:
        file_path: Path to the data file
        
    Returns:
        DataFrame with columns ['angle', 'intensity'], sorted by angle
        
    Raises:
        ValueError: If file cannot be parsed
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Find first line with two numeric values
        data_start = 0
        for i, line in enumerate(lines):
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    float(parts[0])
                    float(parts[1])
                    data_start = i
                    break
                except ValueError:
                    continue
        
        # Read data
        df = pd.read_csv(
            file_path,
            sep=r'\s+',
            skiprows=data_start,
            names=['angle', 'intensity'],
            usecols=[0, 1],
            engine='python'
        )
        
        if len(df) < 10:
            raise ValueError(f"Insufficient data points ({len(df)})")
        
        # Sort by angle
        df = df.sort_values('angle').reset_index(drop=True)
        
        return df
        
    except Exception as e:
        raise ValueError(f"Failed to load {file_path.name}: {e}")


def extract_voltage_from_filename(filename: str) -> Optional[float]:
    """
    Extract voltage in volts from filename (e.g., 'step2_9.5v.txt' -> 9.5).
    
    Args:
        filename: Filename string
        
    Returns:
        Voltage in volts, or None if not found
    """
    # Look for pattern like "9.5v" or "9v"
    match = re.search(r'(\d+(?:\.\d+)?)\s*v(?:_|\.)', filename.lower())
    if match:
        return float(match.group(1))
    return None


def load_vi_table_or_stub(step2_dir: Path) -> pd.DataFrame:
    """
    Load voltage/current table from CSV, or create stub if missing.
    
    Expected file: step2_dir / step2_vi.csv with columns: voltage_V, current_A
    
    Args:
        step2_dir: Path to Step 2 data directory
        
    Returns:
        DataFrame with voltage_V and current_A columns
        
    Raises:
        FileNotFoundError: If stub is created (user must fill in currents)
    """
    vi_csv = step2_dir / "step2_vi.csv"
    
    if vi_csv.exists():
        df = pd.read_csv(vi_csv)
        if 'voltage_V' in df.columns and 'current_A' in df.columns:
            return df
        else:
            raise ValueError(f"CSV {vi_csv} missing required columns (voltage_V, current_A)")
    
    # Create stub: discover voltages from filenames
    step2_files = sorted(step2_dir.glob("step2*.txt"))
    voltages = []
    
    for f in step2_files:
        v = extract_voltage_from_filename(f.name)
        if v is not None:
            voltages.append(v)
    
    if not voltages:
        raise ValueError("Could not extract any voltages from step2 filenames")
    
    # Remove duplicates, sort
    voltages = sorted(set(voltages))
    
    # Create stub CSV
    stub_df = pd.DataFrame({
        'voltage_V': voltages,
        'current_A': [np.nan] * len(voltages)
    })
    
    stub_df.to_csv(vi_csv, index=False)
    
    raise FileNotFoundError(
        f"\nVoltage/current table not found: {vi_csv}\n"
        f"Created stub with {len(voltages)} voltages.\n"
        f"Please fill in the current_A column from LabVIEW and re-run.\n"
        f"Columns: voltage_V, current_A"
    )


# ============================================================================
# UNCERTAINTY PROPAGATION
# ============================================================================

def propagate_theta_peak_uncertainty() -> float:
    """
    Combine instrument angle resolution and peak-picking uncertainty.
    
    Assumes independent contributions add in quadrature:
        u(theta) = sqrt(u_instrument^2 + u_peak_picking^2)
    
    Returns:
        Total uncertainty in degrees
    """
    u_theta = np.sqrt(
        Config.U_THETA_INSTRUMENT**2 + 
        Config.U_THETA_PEAK_PICKING**2
    )
    return float(u_theta)


def lambda_derivative(delta_theta_deg: float, h: float = 1e-4) -> float:
    """
    Compute d(lambda)/d(theta) using finite difference.
    
    Args:
        delta_theta_deg: Angular separation (theta_init - theta_peak) in degrees
        h: Step size for finite difference
        
    Returns:
        Derivative d(lambda)/d(theta) in nm/degree
    """
    try:
        lambda_plus = lambda_from_theta(delta_theta_deg + h)
        lambda_minus = lambda_from_theta(delta_theta_deg - h)
        derivative = (lambda_plus - lambda_minus) / (2 * h)
        return float(derivative)
    except Exception:
        # Fallback: approximate using analytic derivative of prism equation
        # lambda ≈ 10-50 nm in this range; use average sensitivity ~1 nm/deg
        return 2.0  # Conservative estimate


def propagate_lambda_uncertainty(
    delta_theta_deg: float,
    u_theta_deg: float
) -> float:
    """
    Propagate angle uncertainty to wavelength uncertainty.
    
    u(lambda) = |d(lambda)/d(theta)| * u(theta)
    
    Args:
        delta_theta_deg: Angular separation in degrees
        u_theta_deg: Angle uncertainty in degrees
        
    Returns:
        Wavelength uncertainty in nm
    """
    d_lambda_d_theta = lambda_derivative(delta_theta_deg)
    u_lambda = abs(d_lambda_d_theta) * u_theta_deg
    return float(u_lambda)


def propagate_T_uncertainty(
    voltage_V: float,
    current_A: float,
    u_voltage_V: float = Config.U_VOLTAGE,
    u_current_A: float = Config.U_CURRENT,
    T0: float = Config.T0,
    u_T0: float = Config.U_T0,
    R0: float = Config.R0,
    u_R0: float = Config.U_R0,
    alpha0: float = Config.ALPHA0,
    u_alpha0: float = Config.U_ALPHA0
) -> float:
    """
    Propagate V, I, and reference uncertainties to T uncertainty.
    
    Temperature model: T = T0 + ((V/I)/R0 - 1) / alpha0
    
    Partial derivatives (using chain rule):
        ∂T/∂V = 1 / (R0 * alpha0 * I)
        ∂T/∂I = -V / (R0 * alpha0 * I²)
        ∂T/∂T0 = 1
        ∂T/∂R0 = -((V/I) / (R0² * alpha0))
        ∂T/∂alpha0 = -((V/I)/R0 - 1) / alpha0²
    
    u(T) = sqrt( (∂T/∂V*u_V)² + (∂T/∂I*u_I)² + (∂T/∂T0*u_T0)² + 
                 (∂T/∂R0*u_R0)² + (∂T/∂alpha0*u_alpha0)² )
    
    Args:
        voltage_V: Applied voltage (V)
        current_A: Filament current (A)
        u_voltage_V: Voltage uncertainty (V)
        u_current_A: Current uncertainty (A)
        T0: Reference temperature (K)
        u_T0: Reference temperature uncertainty (K)
        R0: Reference resistance (ohm)
        u_R0: Reference resistance uncertainty (ohm)
        alpha0: Temperature coefficient (1/K)
        u_alpha0: Temperature coefficient uncertainty (1/K)
        
    Returns:
        Temperature uncertainty in K
    """
    # Compute partial derivatives
    dT_dV = 1.0 / (R0 * alpha0 * current_A)
    dT_dI = -voltage_V / (R0 * alpha0 * current_A**2)
    dT_dT0 = 1.0
    dT_dR0 = -((voltage_V / current_A) / (R0**2 * alpha0))
    dT_dalpha0 = -((voltage_V / current_A) / R0 - 1.0) / (alpha0**2)
    
    # Combine in quadrature
    u_T_squared = (
        (dT_dV * u_voltage_V)**2 +
        (dT_dI * u_current_A)**2 +
        (dT_dT0 * u_T0)**2 +
        (dT_dR0 * u_R0)**2 +
        (dT_dalpha0 * u_alpha0)**2
    )
    
    u_T = np.sqrt(u_T_squared)
    return float(u_T)


def propagate_inv_T_uncertainty(temperature_K: float, u_T: float) -> float:
    """
    Propagate temperature uncertainty to 1/T uncertainty.
    
    u(1/T) = u(T) / T²
    
    Args:
        temperature_K: Temperature in K
        u_T: Temperature uncertainty in K
        
    Returns:
        Uncertainty in 1/T in 1/K
    """
    u_inv_T = u_T / (temperature_K**2)
    return float(u_inv_T)


def propagate_b_uncertainty(
    lambda_max_nm: float,
    u_lambda_nm: float,
    temperature_K: float,
    u_T_K: float
) -> float:
    """
    Propagate lambda and T uncertainties to Wien constant b uncertainty.
    
    b = lambda_max * T (with lambda converted to m)
    u(b) = sqrt( (T * u_lambda_m)² + (lambda_m * u_T)² )
    
    Args:
        lambda_max_nm: Peak wavelength in nm
        u_lambda_nm: Wavelength uncertainty in nm
        temperature_K: Temperature in K
        u_T_K: Temperature uncertainty in K
        
    Returns:
        Wien constant uncertainty in m·K
    """
    lambda_m = lambda_max_nm * 1e-9
    u_lambda_m = u_lambda_nm * 1e-9
    
    u_b_squared = (temperature_K * u_lambda_m)**2 + (lambda_m * u_T_K)**2
    u_b = np.sqrt(u_b_squared)
    return float(u_b)


# ============================================================================
# PEAK DETECTION
# ============================================================================


def find_theta_peak(
    df: pd.DataFrame,
    angle_window: Tuple[float, float] = Config.PEAK_ANGLE_WINDOW,
    prominence_threshold: float = Config.PROMINENCE_THRESHOLD
) -> float:
    """
    Find the main spectrum peak angle in the specified window.
    
    Args:
        df: DataFrame with 'angle' and 'intensity' columns
        angle_window: (min_angle, max_angle) in degrees to search
        prominence_threshold: Minimum peak prominence
        
    Returns:
        Angle in degrees of the detected peak
        
    Raises:
        ValueError: If peak cannot be found
    """
    angles = df['angle'].values
    intensities = df['intensity'].values
    
    # Filter to window
    mask = (angles >= angle_window[0]) & (angles <= angle_window[1])
    window_angles = angles[mask]
    window_intensities = intensities[mask]
    
    if len(window_angles) < 5:
        raise ValueError(f"Insufficient data points in angle window {angle_window}")
    
    # Find max in window (simple approach)
    max_idx = np.argmax(window_intensities)
    peak_angle = float(window_angles[max_idx])
    
    # Optional: refine with quadratic fit around max
    if len(window_angles) >= 7:
        fit_indices = max(0, max_idx - 3) if max_idx >= 3 else 0
        fit_max = min(len(window_angles), max_idx + 4) if max_idx < len(window_angles) - 4 else len(window_angles)
        
        fit_angles = window_angles[fit_indices:fit_max]
        fit_intensities = window_intensities[fit_indices:fit_max]
        
        # Quadratic fit: y = a*x^2 + b*x + c
        # Vertex at x = -b/(2a)
        try:
            coeffs = np.polyfit(fit_angles, fit_intensities, 2)
            if coeffs[0] != 0:  # Ensure parabola opens downward
                vertex_angle = -coeffs[1] / (2 * coeffs[0])
                # Validate vertex is in reasonable range
                if angle_window[0] <= vertex_angle <= angle_window[1]:
                    peak_angle = vertex_angle
        except Exception:
            pass  # Keep simple max
    
    return peak_angle


# ============================================================================
# WIEN'S LAW CALCULATIONS
# ============================================================================

def lambda_from_theta(delta_theta_deg: float) -> float:
    """
    Convert angular separation to wavelength using Wien's displacement law apparatus equation.
    
    Formula (with theta in radians):
        lambda = sqrt( A / ( sqrt(( (2/sqrt(3))*sin(theta) + 1/2 )^2 + 3/4 ) - B ) )
    
    where A = 13900 nm^2, B = 1.689
    
    Args:
        delta_theta_deg: Angular separation in degrees (theta_init - theta_peak)
        
    Returns:
        Peak wavelength in nm
        
    Raises:
        ValueError: If calculation fails (domain error)
    """
    theta_rad = np.radians(delta_theta_deg)
    
    sin_term = (2.0 / np.sqrt(3.0)) * np.sin(theta_rad) + 0.5
    sqrt_term = np.sqrt(sin_term**2 + 0.75)
    denominator = sqrt_term - Config.WIEN_B
    
    if denominator <= 0:
        raise ValueError(f"Invalid wavelength calculation: denominator = {denominator} (delta_theta = {delta_theta_deg}°)")
    
    lambda_nm = np.sqrt(Config.WIEN_A / denominator)
    
    if not np.isfinite(lambda_nm):
        raise ValueError(f"Non-finite wavelength: {lambda_nm}")
    
    return float(lambda_nm)


def compute_temperature(voltage_V: float, current_A: float) -> float:
    """
    Compute filament temperature from voltage and current.
    
    T = T0 + ( (V/I)/R0 - 1 ) / alpha0
    
    where:
      T0 = 293 K (reference)
      R0 = 1.1 ohm (reference resistance)
      alpha0 = 4.5e-3 1/K (temperature coefficient)
    
    Args:
        voltage_V: Applied voltage (V)
        current_A: Filament current (A)
        
    Returns:
        Temperature in Kelvin
        
    Raises:
        ValueError: If calculation fails
    """
    if current_A <= 0 or voltage_V < 0:
        raise ValueError(f"Invalid V/I: V={voltage_V}, I={current_A}")
    
    resistance = voltage_V / current_A
    T = Config.T0 + ((resistance / Config.R0) - 1.0) / Config.ALPHA0
    
    if T < 0 or not np.isfinite(T):
        raise ValueError(f"Invalid temperature: {T} K")
    
    return float(T)


def compute_wien_constant(lambda_nm: float, temperature_K: float) -> float:
    """
    Compute Wien constant: b = lambda_max * T (in m·K).
    
    Args:
        lambda_nm: Peak wavelength in nm
        temperature_K: Temperature in K
        
    Returns:
        Wien constant in m·K
    """
    lambda_m = lambda_nm * 1e-9
    b = lambda_m * temperature_K
    return float(b)


# ============================================================================
# PLOTTING
# ============================================================================

def plot_lambda_vs_inv_t(results_df: pd.DataFrame, output_dir: Path) -> None:
    """
    Plot lambda_max vs 1/T with error bars, weighted regression, residuals, and chi-squared.
    
    Args:
        results_df: Results DataFrame with columns including u_lambda_max_nm, u_inv_T
        output_dir: Output directory
    """
    # Use columns with uncertainties if available, otherwise fall back
    valid_cols = ['lambda_max_nm', 'T_K', 'u_lambda_max_nm', 'u_inv_T']
    valid_df = results_df.dropna(subset=['lambda_max_nm', 'T_K'])
    
    if valid_df.empty:
        return
    
    inv_t = 1.0 / valid_df['T_K'].values
    lambda_nm = valid_df['lambda_max_nm'].values
    
    # Get uncertainties if available
    if 'u_lambda_max_nm' in valid_df.columns and 'u_inv_T' in valid_df.columns:
        u_lambda_nm = valid_df['u_lambda_max_nm'].fillna(0).values
        u_inv_t = valid_df['u_inv_T'].fillna(0).values
    else:
        u_lambda_nm = np.zeros_like(lambda_nm)
        u_inv_t = np.zeros_like(inv_t)
    
    # Create figure with two panels
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 9), sharex=True, 
                                    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
    
    # Top panel: data and fit with error bars
    ax1.errorbar(inv_t, lambda_nm, xerr=u_inv_t, yerr=u_lambda_nm, 
                fmt='o', markersize=8, capsize=5, alpha=0.6, 
                color='tab:blue', label='Data (±1σ)', zorder=3)
    
    # Weighted linear regression
    chi_squared = np.nan
    reduced_chi_sq = np.nan
    slope_unc = np.nan
    intercept_unc = np.nan
    
    try:
        # Compute weights: w = 1/u_y^2 (weights from y-uncertainties)
        # Handle zero uncertainties with small default
        u_lambda_safe = np.where(u_lambda_nm > 0, u_lambda_nm, np.median(u_lambda_nm[u_lambda_nm > 0]) if np.any(u_lambda_nm > 0) else 1.0)
        weights = 1.0 / (u_lambda_safe**2)
        
        # Weighted polyfit
        coeffs = np.polyfit(inv_t, lambda_nm, 1, w=weights)
        fit_line = np.polyval(coeffs, inv_t)
        residuals = lambda_nm - fit_line
        
        # Compute chi-squared with weights
        chi_squared = np.sum(weights * (residuals**2))
        dof = len(lambda_nm) - 2  # 2 parameters in linear fit
        reduced_chi_sq = chi_squared / dof
        
        # Estimate parameter uncertainties from covariance matrix
        # For weighted fit: Cov = (X^T W X)^-1
        n = len(inv_t)
        X = np.vstack([inv_t, np.ones(n)]).T
        W = np.diag(weights)
        XtWX = X.T @ W @ X
        try:
            cov_matrix = np.linalg.inv(XtWX)
            slope_unc = np.sqrt(cov_matrix[0, 0])
            intercept_unc = np.sqrt(cov_matrix[1, 1])
        except np.linalg.LinAlgError:
            slope_unc = np.nan
            intercept_unc = np.nan
        
        slope = coeffs[0]
        intercept = coeffs[1]
        
        # Plot fit line
        if not np.isnan(slope_unc) and not np.isnan(intercept_unc):
            label_str = (f'Weighted linear fit\n'
                        f'Slope: ({slope:.2f} ± {slope_unc:.2f}) nm·K\n'
                        f'Intercept: ({intercept:.1f} ± {intercept_unc:.1f}) nm')
        else:
            label_str = f'Linear fit: y={slope:.1f}x+{intercept:.1f}'
        
        ax1.plot(inv_t, fit_line, 'r--', linewidth=2, label=label_str)
        
        # Display chi-squared statistics in red box
        ax1.text(0.98, 0.95, f'χ² = {chi_squared:.2f}\nχ²_red = {reduced_chi_sq:.2f}', 
                transform=ax1.transAxes, va='top', ha='right', fontsize=11, color='red', 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='red', linewidth=1.5))
        
        # Bottom panel: residuals with error bars
        ax2.errorbar(inv_t, residuals, yerr=u_lambda_nm, fmt='o', markersize=6, 
                    capsize=4, alpha=0.6, color='tab:blue', zorder=3)
        ax2.axhline(0, color='red', linestyle='--', linewidth=1.5, zorder=2)
        
        # Shaded band for ±1σ of residuals
        sigma_resid = np.std(residuals, ddof=1)
        ax2.fill_between([inv_t.min(), inv_t.max()], -sigma_resid, sigma_resid, 
                        color='red', alpha=0.15, label=f'±1σ residuals ({sigma_resid:.1f} nm)', zorder=1)
        
        ax2.set_ylabel('Residuals (nm)', fontsize=10)
        ax2.grid(True, alpha=0.3)
        ax2.legend(loc='best', fontsize=9)
        
    except Exception as e:
        ax1.text(0.05, 0.95, f'Fit failed: {e}', transform=ax1.transAxes, 
                va='top', fontsize=10, color='red')
        print(f"  Warning: Lambda vs 1/T fit failed: {e}")
    
    ax1.set_ylabel('λ_max (nm)', fontsize=11)
    ax1.set_title('Wien\'s Law: λ_max vs 1/T with Uncertainty Error Bars', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left', fontsize=10)
    
    ax2.set_xlabel('1/T (K⁻¹)', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'step2_lambda_vs_inv_t.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_wien_constant_vs_t(results_df: pd.DataFrame, output_dir: Path) -> None:
    """
    Plot Wien constant vs T to check constancy, with error bars.
    
    Args:
        results_df: Results DataFrame with columns including u_b_mK
        output_dir: Output directory
    """
    valid_df = results_df.dropna(subset=['b_mK', 'T_K'])
    if valid_df.empty:
        return
    
    t_k = valid_df['T_K'].values
    b_mk = valid_df['b_mK'].values
    
    # Get uncertainties if available
    if 'u_b_mK' in valid_df.columns:
        u_b_mk = valid_df['u_b_mK'].fillna(0).values
    else:
        u_b_mk = np.zeros_like(b_mk)
    
    fig, ax = plt.subplots(figsize=(11, 7))
    
    # Plot data with error bars
    ax.errorbar(t_k, b_mk, yerr=u_b_mk, fmt='o', markersize=9, capsize=6, 
               alpha=0.7, color='tab:green', label='Measurements (±1σ)', zorder=3)
    
    # Literature value with horizontal line
    ax.axhline(Config.WIEN_B0, color='red', linestyle='--', linewidth=2, 
              label=f'Literature b₀ = {Config.WIEN_B0:.4f} m·K', zorder=2)
    
    # Compute and show mean ± std of measurements
    mean_b = np.mean(b_mk)
    std_b = np.std(b_mk, ddof=1)
    ax.axhline(mean_b, color='orange', linestyle='-.', linewidth=1.5, 
              label=f'Measurement mean = {mean_b:.4e} m·K', zorder=2, alpha=0.7)
    ax.fill_between([t_k.min(), t_k.max()], mean_b - std_b, mean_b + std_b, 
                    color='orange', alpha=0.15, zorder=1, label=f'±1σ range = ±{std_b:.4e} m·K')
    
    ax.set_xlabel('T (K)', fontsize=11)
    ax.set_ylabel('b (m·K)', fontsize=11)
    ax.set_title('Wien Constant vs Temperature (Testing Constancy)', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc='best')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    
    plt.tight_layout()
    plt.savefig(output_dir / 'step2_wien_constant_vs_t.png', dpi=300, bbox_inches='tight')
    plt.close()



# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def analyze_step2_file(
    file_path: Path,
    vi_df: pd.DataFrame,
    config: Config,
    output_dir: Path
) -> Dict:
    """
    Analyze a single Step 2 file with full uncertainty propagation.
    
    Args:
        file_path: Path to Step 2 data file
        vi_df: Voltage/current DataFrame
        config: Configuration object
        output_dir: Output directory
        
    Returns:
        Dictionary with analysis results and uncertainties
    """
    result = {
        'file': file_path.name,
        'voltage_V': np.nan,
        'current_A': np.nan,
        'theta_peak_deg': np.nan,
        'u_theta_peak_deg': np.nan,
        'delta_theta_deg': np.nan,
        'lambda_max_nm': np.nan,
        'u_lambda_max_nm': np.nan,
        'T_K': np.nan,
        'u_T_K': np.nan,
        'inv_T': np.nan,
        'u_inv_T': np.nan,
        'b_mK': np.nan,
        'u_b_mK': np.nan,
        'notes': ''
    }
    
    try:
        # Extract voltage from filename
        voltage = extract_voltage_from_filename(file_path.name)
        if voltage is None:
            raise ValueError("Could not extract voltage from filename")
        result['voltage_V'] = voltage
        
        # Look up current
        vi_row = vi_df[np.isclose(vi_df['voltage_V'], voltage, atol=0.01)]
        if vi_row.empty:
            raise ValueError(f"Voltage {voltage}V not found in VI table")
        
        current = vi_row.iloc[0]['current_A']
        if pd.isna(current):
            raise ValueError(f"Current for {voltage}V is NaN")
        result['current_A'] = current
        
        # Load and analyze scan
        df = load_scan(file_path)
        theta_peak = find_theta_peak(df, angle_window=config.PEAK_ANGLE_WINDOW)
        result['theta_peak_deg'] = theta_peak
        
        # Compute angle uncertainty
        u_theta_peak = propagate_theta_peak_uncertainty()
        result['u_theta_peak_deg'] = u_theta_peak
        
        # Compute delta_theta
        delta_theta = config.THETA_INIT_DEG - theta_peak
        result['delta_theta_deg'] = delta_theta
        
        # Compute wavelength with uncertainty
        lambda_nm = lambda_from_theta(delta_theta)
        result['lambda_max_nm'] = lambda_nm
        
        u_lambda_nm = propagate_lambda_uncertainty(delta_theta, u_theta_peak)
        result['u_lambda_max_nm'] = u_lambda_nm
        
        # Compute temperature with uncertainty
        temp_K = compute_temperature(voltage, current)
        result['T_K'] = temp_K
        
        u_temp_K = propagate_T_uncertainty(voltage, current)
        result['u_T_K'] = u_temp_K
        
        # Compute 1/T and its uncertainty
        inv_T = 1.0 / temp_K
        result['inv_T'] = inv_T
        
        u_inv_T = propagate_inv_T_uncertainty(temp_K, u_temp_K)
        result['u_inv_T'] = u_inv_T
        
        # Compute Wien constant with uncertainty
        b_mk = compute_wien_constant(lambda_nm, temp_K)
        result['b_mK'] = b_mk
        
        u_b_mk = propagate_b_uncertainty(lambda_nm, u_lambda_nm, temp_K, u_temp_K)
        result['u_b_mK'] = u_b_mk
        
        result['notes'] = 'OK'
        
    except Exception as e:
        result['notes'] = f"ERROR: {e}"
        warnings.warn(f"Failed to analyze {file_path.name}: {e}")
    
    return result


def main():
    """Main entry point for Step 2 Wien analysis."""
    parser = argparse.ArgumentParser(
        description='Analyze Step 2 data to extract Wien constant',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--data_dir', type=str, default=Config.DATA_DIR,
                        help='Base data directory')
    parser.add_argument('--out_dir', type=str, default=Config.OUTPUT_DIR,
                        help='Output directory')
    parser.add_argument('--theta_init', type=float, default=Config.THETA_INIT_DEG,
                        help='Theta_init from Step 1 (degrees)')
    parser.add_argument('--debug', action='store_true',
                        help='Enable debug mode')
    
    args = parser.parse_args()
    
    # Update config
    config = Config()
    config.THETA_INIT_DEG = args.theta_init
    config.DEBUG = args.debug
    
    step2_dir = Path(args.data_dir) / "step2"
    if not step2_dir.exists():
        step2_dir = Path(args.data_dir)
    if not step2_dir.exists():
        print(f"ERROR: Step 2 directory not found: {step2_dir}")
        return
    
    # Load VI table
    try:
        vi_df = load_vi_table_or_stub(step2_dir)
        print(f"Loaded VI table with {len(vi_df)} voltage points")
    except FileNotFoundError as e:
        print(str(e))
        return
    except Exception as e:
        print(f"ERROR loading VI table: {e}")
        return
    
    # Find Step 2 files
    step2_files = sorted(step2_dir.glob("step2*.txt"))
    if not step2_files:
        print(f"ERROR: No Step 2 files found in {step2_dir}")
        return
    
    print(f"Found {len(step2_files)} Step 2 files")
    print("=" * 70)
    
    # Create output directory
    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Analyze each file
    results = []
    for file_path in step2_files:
        print(f"Processing: {file_path.name}")
        result = analyze_step2_file(file_path, vi_df, config, output_dir)
        results.append(result)
        
        if not np.isnan(result['T_K']):
            print(f"  → V={result['voltage_V']:.1f}V, I={result['current_A']:.3f}A, "
                  f"T={result['T_K']:.0f}K, λ={result['lambda_max_nm']:.1f}nm, "
                  f"b={result['b_mK']*1e3:.3f}×10⁻³ m·K")
        else:
            print(f"  → {result['notes']}")
    
    # Save results table
    results_df = pd.DataFrame(results)
    results_path = output_dir / 'step2_wien_results.csv'
    results_df.to_csv(results_path, index=False)
    print(f"\nResults saved to: {results_path}")
    
    # Generate plots
    plot_lambda_vs_inv_t(results_df, output_dir)
    plot_wien_constant_vs_t(results_df, output_dir)
    print(f"Plots saved to: {output_dir}/")
    
    # Compute statistics
    valid_b = results_df['b_mK'].dropna().values
    if len(valid_b) > 0:
        mean_b = np.mean(valid_b)
        std_b = np.std(valid_b, ddof=1) if len(valid_b) > 1 else 0.0
        percent_diff = ((mean_b - Config.WIEN_B0) / Config.WIEN_B0) * 100
        
        summary_lines = []
        summary_lines.append("=" * 70)
        summary_lines.append("STEP 2 ANALYSIS SUMMARY: Wien's Displacement Law")
        summary_lines.append("=" * 70)
        summary_lines.append(f"Files processed: {len(results)}")
        summary_lines.append(f"Successful analyses: {len(valid_b)}")
        summary_lines.append("")
        summary_lines.append(f"Wien constant mean:              {mean_b*1e3:.4f}×10⁻³ m·K")
        summary_lines.append(f"Wien constant std dev:           {std_b*1e3:.4f}×10⁻³ m·K")
        summary_lines.append(f"Literature b₀:                   {Config.WIEN_B0*1e3:.4f}×10⁻³ m·K")
        summary_lines.append(f"Percent difference:              {percent_diff:.2f}%")
        summary_lines.append("=" * 70)
        
        summary_text = "\n".join(summary_lines)
        
        # Save and print
        summary_path = output_dir / 'step2_wien_summary.txt'
        with open(summary_path, 'w') as f:
            f.write(summary_text)
        
        print("\n" + summary_text)
        print(f"\nSummary saved to: {summary_path}")
    
    print("\n✓ Step 2 analysis complete!")


if __name__ == '__main__':
    main()
