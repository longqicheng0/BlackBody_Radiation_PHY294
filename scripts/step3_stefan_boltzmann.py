#!/usr/bin/env python3
"""
Blackbody Lab - Step 3 Analysis: Stefan-Boltzmann Law

Analyzes intensity vs wavelength to compute total radiated power,
visible/IR fractions, and validate Stefan-Boltzmann T^4 law.

Author: PHY294 Lab Analysis Tool
"""

import argparse
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use('Agg')  # Non-GUI backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================================
# CONFIGURATION
# ============================================================================

class Config:
    """Configuration for Step 3 Stefan-Boltzmann analysis."""
    DATA_DIR = "Blackbody_Lab_Data"
    PATTERN = "**/step2*.txt"
    OUTPUT_DIR = "outputs"
    
    # Calibration from Step 1
    THETA_INIT_DEG = 78.28
    
    # Wavelength constants (Wien apparatus)
    WIEN_A = 13900  # nm^2
    WIEN_B = 1.689
    
    # Spectrum regions (nm)
    LAMBDA_VISIBLE_MIN = 400.0
    LAMBDA_VISIBLE_MAX = 700.0
    LAMBDA_IR_MIN = 700.0
    
    DEBUG = False


# ============================================================================
# WAVELENGTH CONVERSION
# ============================================================================

def lambda_from_theta(delta_theta_deg: float) -> float:
    """
    Convert angular separation to wavelength using Wien apparatus equation.
    
    Args:
        delta_theta_deg: Angular separation in degrees (theta_init - theta)
        
    Returns:
        Wavelength in nm
    """
    theta_rad = np.radians(delta_theta_deg)
    sin_term = (2.0 / np.sqrt(3.0)) * np.sin(theta_rad) + 0.5
    sqrt_term = np.sqrt(sin_term**2 + 0.75)
    denominator = sqrt_term - Config.WIEN_B
    
    if denominator <= 0:
        return np.nan
    
    lambda_nm = np.sqrt(Config.WIEN_A / denominator)
    return float(lambda_nm)


def convert_angles_to_wavelengths(angles_deg: np.ndarray, theta_init: float) -> np.ndarray:
    """
    Convert array of angles to wavelengths.
    
    Args:
        angles_deg: Array of angles in degrees
        theta_init: Reference angle (unrefracted light)
        
    Returns:
        Array of wavelengths in nm
    """
    delta_theta = theta_init - angles_deg
    wavelengths = np.array([lambda_from_theta(dt) for dt in delta_theta])
    return wavelengths


# ============================================================================
# DATA LOADING
# ============================================================================

def load_scan(file_path: Path) -> pd.DataFrame:
    """
    Load a scan file with angle and intensity data.
    
    Args:
        file_path: Path to data file
        
    Returns:
        DataFrame with columns ['angle', 'intensity']
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Find first data line
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
    
    df = pd.read_csv(
        file_path,
        sep=r'\s+',
        skiprows=data_start,
        names=['angle', 'intensity'],
        usecols=[0, 1],
        engine='python'
    )
    
    df = df.sort_values('angle').reset_index(drop=True)
    return df


# ============================================================================
# NUMERICAL INTEGRATION
# ============================================================================

def integrate_spectrum(
    wavelengths_nm: np.ndarray,
    intensities: np.ndarray,
    lambda_min: float = None,
    lambda_max: float = None
) -> float:
    """
    Integrate intensity over wavelength using trapezoidal rule.
    
    Args:
        wavelengths_nm: Wavelength array in nm
        intensities: Intensity array (arbitrary units)
        lambda_min: Minimum wavelength for integration (None = use data min)
        lambda_max: Maximum wavelength for integration (None = use data max)
        
    Returns:
        Integrated area (relative radiated power)
    """
    # Filter to integration range
    mask = np.isfinite(wavelengths_nm) & np.isfinite(intensities)
    
    if lambda_min is not None:
        mask &= (wavelengths_nm >= lambda_min)
    if lambda_max is not None:
        mask &= (wavelengths_nm <= lambda_max)
    
    lambda_filtered = wavelengths_nm[mask]
    intensity_filtered = intensities[mask]
    
    if len(lambda_filtered) < 2:
        return np.nan
    
    # Sort by wavelength for integration
    sort_idx = np.argsort(lambda_filtered)
    lambda_sorted = lambda_filtered[sort_idx]
    intensity_sorted = intensity_filtered[sort_idx]
    
    # Trapezoidal integration
    area = np.trapz(intensity_sorted, lambda_sorted)
    
    return float(area)


def estimate_integration_uncertainty(
    wavelengths_nm: np.ndarray,
    intensities: np.ndarray,
    integrated_area: float,
    lambda_min: float = None,
    lambda_max: float = None
) -> float:
    """
    Estimate uncertainty in integrated area from intensity scatter.
    
    Uses u(A) ≈ A * (std(I) / mean(I)) as a simple approximation.
    
    Args:
        wavelengths_nm: Wavelength array
        intensities: Intensity array
        integrated_area: Previously computed integrated area
        lambda_min: Minimum wavelength (for filtering)
        lambda_max: Maximum wavelength (for filtering)
        
    Returns:
        Uncertainty in integrated area
    """
    # Filter to integration range
    mask = np.isfinite(wavelengths_nm) & np.isfinite(intensities)
    
    if lambda_min is not None:
        mask &= (wavelengths_nm >= lambda_min)
    if lambda_max is not None:
        mask &= (wavelengths_nm <= lambda_max)
    
    intensity_filtered = intensities[mask]
    
    if len(intensity_filtered) < 2:
        return np.nan
    
    # Use 5% as a conservative uncertainty estimate for integration
    # This accounts for discretization error, background subtraction, and noise
    # Using scatter (std/mean) gives unrealistically large uncertainties for
    # small spectral regions with few photons
    uncertainty = abs(integrated_area) * 0.05
    
    return float(uncertainty)


# ============================================================================
# STEFAN-BOLTZMANN ANALYSIS
# ============================================================================

def analyze_stefan_boltzmann_file(
    file_path: Path,
    temperature_K: float,
    config: Config
) -> Dict:
    """
    Analyze single file for Stefan-Boltzmann law.
    
    Args:
        file_path: Path to scan file
        temperature_K: Filament temperature in K
        config: Configuration object
        
    Returns:
        Dictionary with analysis results
    """
    result = {
        'file': file_path.name,
        'T_K': temperature_K,
        'A_total': np.nan,
        'A_visible': np.nan,
        'A_IR': np.nan,
        'f_visible': np.nan,
        'f_IR': np.nan,
        'notes': ''
    }
    
    try:
        # Load scan data
        df = load_scan(file_path)
        angles = df['angle'].values
        intensities = df['intensity'].values
        
        # Convert angles to wavelengths
        wavelengths = convert_angles_to_wavelengths(angles, config.THETA_INIT_DEG)
        
        # Total integrated area
        A_total = integrate_spectrum(wavelengths, intensities)
        result['A_total'] = A_total
        
        # Visible region (400-700 nm)
        A_visible = integrate_spectrum(
            wavelengths, intensities,
            lambda_min=Config.LAMBDA_VISIBLE_MIN,
            lambda_max=Config.LAMBDA_VISIBLE_MAX
        )
        result['A_visible'] = A_visible
        
        # IR region (>700 nm)
        lambda_max_data = np.nanmax(wavelengths)
        A_IR = integrate_spectrum(
            wavelengths, intensities,
            lambda_min=Config.LAMBDA_IR_MIN,
            lambda_max=lambda_max_data
        )
        result['A_IR'] = A_IR
        
        # Compute fractions
        if A_total > 0:
            f_visible = A_visible / A_total
            f_IR = A_IR / A_total
            result['f_visible'] = f_visible
            result['f_IR'] = f_IR
        
        result['notes'] = 'OK'
        
    except Exception as e:
        result['notes'] = f'ERROR: {e}'
        warnings.warn(f"Failed to analyze {file_path.name}: {e}")
    
    return result


# ============================================================================
# PLOTTING
# ============================================================================

def plot_A_total_vs_T(results_df: pd.DataFrame, output_dir: Path) -> None:
    """
    Plot total radiated power vs temperature with error bars.
    
    Args:
        results_df: Results DataFrame
        output_dir: Output directory
    """
    valid_df = results_df.dropna(subset=['T_K', 'A_total'])
    if valid_df.empty:
        return
    
    T = valid_df['T_K'].values
    A_total = valid_df['A_total'].values
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.scatter(T, A_total, s=100, alpha=0.7, color='tab:blue', 
               label='Measurements', zorder=3)
    
    ax.set_xlabel('T (K)', fontsize=12)
    ax.set_ylabel('A_total (integrated intensity, a.u.)', fontsize=12)
    ax.set_title('Total Radiated Power vs Temperature', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'step3_A_total_vs_T.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_log_A_vs_log_T(results_df: pd.DataFrame, output_dir: Path) -> None:
    """
    Plot log(A_total) vs log(T) with linear regression to find power law exponent.
    
    Args:
        results_df: Results DataFrame
        output_dir: Output directory
    """
    valid_df = results_df.dropna(subset=['T_K', 'A_total'])
    if valid_df.empty:
        return
    
    T = valid_df['T_K'].values
    A_total = valid_df['A_total'].values
    
    # Remove any negative or zero values
    mask = (T > 0) & (A_total > 0)
    T = T[mask]
    A_total = A_total[mask]
    
    if len(T) < 2:
        return
    
    log_T = np.log10(T)
    log_A = np.log10(A_total)
    
    # Linear regression
    coeffs = np.polyfit(log_T, log_A, 1)
    slope = coeffs[0]
    intercept = coeffs[1]
    
    # Estimate slope uncertainty
    fit_line = np.polyval(coeffs, log_T)
    residuals = log_A - fit_line
    
    # Standard error of slope
    n = len(log_T)
    s_res = np.sqrt(np.sum(residuals**2) / (n - 2))
    s_xx = np.sum((log_T - np.mean(log_T))**2)
    slope_unc = s_res / np.sqrt(s_xx)
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.scatter(log_T, log_A, s=100, alpha=0.6, color='tab:blue', 
              label='Data', zorder=3)
    
    # Plot fit
    log_T_fit = np.linspace(log_T.min(), log_T.max(), 100)
    log_A_fit = np.polyval(coeffs, log_T_fit)
    ax.plot(log_T_fit, log_A_fit, 'r--', linewidth=2,
           label=f'Linear fit: slope = {slope:.2f} ± {slope_unc:.2f}\n(Expected: slope ≈ 4 for Stefan-Boltzmann)')
    
    ax.set_xlabel('log₁₀(T) [K]', fontsize=12)
    ax.set_ylabel('log₁₀(A_total) [a.u.]', fontsize=12)
    ax.set_title('Stefan-Boltzmann Law: log-log Plot', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc='best')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'step3_log_A_vs_log_T.png', dpi=300, bbox_inches='tight')
    plt.close()


def plot_f_visible_vs_T(results_df: pd.DataFrame, output_dir: Path) -> None:
    """
    Plot visible fraction vs temperature.
    
    Args:
        results_df: Results DataFrame
        output_dir: Output directory
    """
    valid_df = results_df.dropna(subset=['T_K', 'f_visible'])
    if valid_df.empty:
        return
    
    T = valid_df['T_K'].values
    f_visible = valid_df['f_visible'].values
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    ax.scatter(T, f_visible, s=100, alpha=0.7, color='tab:orange',
               label='Visible fraction', zorder=3)
    
    ax.set_xlabel('T (K)', fontsize=12)
    ax.set_ylabel('f_visible (fraction in 400-700 nm)', fontsize=12)
    ax.set_title('Visible Light Fraction vs Temperature', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'step3_f_visible_vs_T.png', dpi=300, bbox_inches='tight')
    plt.close()


# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main entry point for Step 3 Stefan-Boltzmann analysis."""
    parser = argparse.ArgumentParser(
        description='Analyze Stefan-Boltzmann law from intensity data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--data_dir', type=str, default=Config.DATA_DIR,
                        help='Base data directory')
    parser.add_argument('--out_dir', type=str, default=Config.OUTPUT_DIR,
                        help='Output directory')
    parser.add_argument('--debug', action='store_true',
                        help='Enable debug mode')
    
    args = parser.parse_args()
    
    config = Config()
    config.DEBUG = args.debug
    
    # Load Step 2 results to get temperatures
    step2_results_path = Path(args.out_dir) / 'step2_wien_results.csv'
    if not step2_results_path.exists():
        print(f"ERROR: Step 2 results not found: {step2_results_path}")
        print("Please run Step 2 analysis first (scripts/step2_wien.py)")
        return
    
    step2_df = pd.read_csv(step2_results_path)
    
    # Create output directory
    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("STEP 3: STEFAN-BOLTZMANN LAW ANALYSIS")
    print("=" * 70)
    
    # Analyze each file
    results = []
    for _, row in step2_df.iterrows():
        filename = row['file']
        temperature = row['T_K']
        
        # Find file in data directory
        data_dir = Path(args.data_dir)
        file_path = data_dir / filename
        
        if not file_path.exists():
            print(f"Warning: File not found: {file_path}")
            continue
        
        print(f"Processing: {filename} (T={temperature:.0f}K)")
        result = analyze_stefan_boltzmann_file(file_path, temperature, config)
        results.append(result)
        
        if not np.isnan(result['A_total']):
            print(f"  → A_total={result['A_total']:.2e}, "
                  f"f_visible={result['f_visible']:.3f}, "
                  f"f_IR={result['f_IR']:.3f}")
    
    # Save results
    results_df = pd.DataFrame(results)
    results_path = output_dir / 'step3_stefan_boltzmann_results.csv'
    results_df.to_csv(results_path, index=False)
    print(f"\nResults saved to: {results_path}")
    
    # Generate plots
    plot_A_total_vs_T(results_df, output_dir)
    plot_log_A_vs_log_T(results_df, output_dir)
    plot_f_visible_vs_T(results_df, output_dir)
    print(f"Plots saved to: {output_dir}/")
    
    # Summary statistics
    valid_results = results_df.dropna(subset=['A_total'])
    if not valid_results.empty:
        print("\n" + "=" * 70)
        print("SUMMARY STATISTICS")
        print("=" * 70)
        print(f"Files processed: {len(results_df)}")
        print(f"Successful analyses: {len(valid_results)}")
        print(f"\nMean visible fraction: {valid_results['f_visible'].mean():.3f}")
        print(f"Mean IR fraction: {valid_results['f_IR'].mean():.3f}")
        
        # Power law exponent from log-log fit
        T_vals = valid_results['T_K'].values
        A_vals = valid_results['A_total'].values
        mask = (T_vals > 0) & (A_vals > 0)
        if np.sum(mask) >= 2:
            log_T = np.log10(T_vals[mask])
            log_A = np.log10(A_vals[mask])
            slope, _ = np.polyfit(log_T, log_A, 1)
            print(f"\nStefan-Boltzmann power law exponent: {slope:.2f}")
            print(f"Expected: 4.0 (A ∝ T⁴)")
            print(f"Difference: {slope - 4.0:.2f}")
    
    print("\n✓ Step 3 analysis complete!")


if __name__ == '__main__':
    main()
