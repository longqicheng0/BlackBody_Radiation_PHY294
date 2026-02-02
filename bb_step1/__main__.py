#!/usr/bin/env python3
"""
Blackbody Lab - Step 1 Analysis: theta_init Extraction

This script analyzes Step 1 data files to find the small "direct/unrefracted light" peak
(theta_init) that occurs when light bypasses the prism. It processes all Step 1 files,
generates plots, and computes statistics.

Author: PHY294 Lab Analysis Tool
"""

import argparse
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use('Agg')  # Non-GUI backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Optional scipy for better smoothing
try:
    from scipy.signal import savgol_filter, find_peaks
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    warnings.warn("scipy not available; using fallback methods")


# ============================================================================
# CONFIGURATION DEFAULTS
# ============================================================================

class Config:
    """Default configuration parameters for analysis."""
    # Data discovery
    DATA_DIR = "Blackbody_Lab_Data"
    PATTERN = "**/step1*.txt"
    OUTPUT_DIR = "outputs"
    
    # Preprocessing
    SMOOTH_WINDOW = 11  # Must be odd for Savitzky-Golay
    SMOOTH_POLY = 3     # Polynomial order for Savitzky-Golay
    BASELINE_PERCENTILE = 5  # Percentile for baseline estimation
    
    # Peak detection
    PROMINENCE_THRESHOLD = 0.01  # Minimum prominence (volts)
    MIN_SEPARATION_DEG = 2.0     # Minimum angle separation between main and small peak
    MIN_HEIGHT_RATIO = 0.05      # Small peak must be at least this fraction of main peak
    
    # Debug mode
    DEBUG = False


# ============================================================================
# DATA LOADING AND PARSING
# ============================================================================

def load_scan(file_path: Path) -> pd.DataFrame:
    """
    Load a blackbody scan file (tab-separated).
    
    Automatically skips non-numeric header lines and extracts two columns:
    angle (degrees) and intensity (volts).
    
    Args:
        file_path: Path to the data file
        
    Returns:
        DataFrame with columns ['angle', 'intensity']
        
    Raises:
        ValueError: If file cannot be parsed or has insufficient data
    """
    try:
        # Read file, skipping lines until we find numeric data
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Find first line with two numeric values (tab or space separated)
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
        
        # Read from data_start
        df = pd.read_csv(
            file_path,
            sep=r'\s+',  # Any whitespace
            skiprows=data_start,
            names=['angle', 'intensity'],
            usecols=[0, 1],
            engine='python'
        )
        
        # Validate
        if len(df) < 10:
            raise ValueError(f"Insufficient data points ({len(df)}) in {file_path.name}")
        
        # Sort by angle
        df = df.sort_values('angle').reset_index(drop=True)
        
        return df
        
    except Exception as e:
        raise ValueError(f"Failed to load {file_path.name}: {e}")


# ============================================================================
# PREPROCESSING
# ============================================================================

def preprocess_scan(
    df: pd.DataFrame,
    smooth_window: int = Config.SMOOTH_WINDOW,
    smooth_poly: int = Config.SMOOTH_POLY,
    baseline_percentile: float = Config.BASELINE_PERCENTILE
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Preprocess scan data: smooth and baseline-correct.
    
    Args:
        df: DataFrame with 'angle' and 'intensity' columns
        smooth_window: Window size for smoothing (must be odd)
        smooth_poly: Polynomial order for Savitzky-Golay (if scipy available)
        baseline_percentile: Percentile for baseline estimation
        
    Returns:
        Tuple of (angles, raw_intensity, smoothed_intensity, baseline_corrected)
    """
    angles = df['angle'].values
    raw_intensity = df['intensity'].values
    
    # Smoothing
    if SCIPY_AVAILABLE and len(angles) > smooth_window:
        # Ensure window is odd and valid
        window = smooth_window if smooth_window % 2 == 1 else smooth_window + 1
        window = min(window, len(angles) - 1)
        if window % 2 == 0:
            window -= 1
        window = max(3, window)
        
        try:
            smoothed = savgol_filter(raw_intensity, window, smooth_poly)
        except Exception:
            # Fallback to moving average
            smoothed = moving_average(raw_intensity, window)
    else:
        # Moving average fallback
        window = min(smooth_window, len(angles))
        smoothed = moving_average(raw_intensity, window)
    
    # Baseline correction (subtract low percentile to handle drift)
    baseline = np.percentile(smoothed, baseline_percentile)
    baseline_corrected = smoothed - baseline
    
    return angles, raw_intensity, smoothed, baseline_corrected


def moving_average(data: np.ndarray, window: int) -> np.ndarray:
    """Simple moving average smoothing."""
    if window < 3:
        return data
    
    # Ensure window is odd
    if window % 2 == 0:
        window += 1
    
    half_window = window // 2
    smoothed = np.copy(data)
    
    for i in range(len(data)):
        start = max(0, i - half_window)
        end = min(len(data), i + half_window + 1)
        smoothed[i] = np.mean(data[start:end])
    
    return smoothed


# ============================================================================
# PEAK DETECTION
# ============================================================================

def detect_peaks(
    angles: np.ndarray,
    signal: np.ndarray,
    prominence_threshold: float = Config.PROMINENCE_THRESHOLD
) -> pd.DataFrame:
    """
    Detect peaks in the signal.
    
    Args:
        angles: Angle values (degrees)
        signal: Intensity signal (baseline-corrected)
        prominence_threshold: Minimum prominence for peak detection
        
    Returns:
        DataFrame with columns: ['index', 'angle', 'height', 'prominence']
        sorted by prominence (descending)
    """
    if SCIPY_AVAILABLE:
        # Use scipy's find_peaks
        peak_indices, properties = find_peaks(
            signal,
            prominence=prominence_threshold,
            height=0  # Only positive peaks after baseline correction
        )
        
        if len(peak_indices) == 0:
            return pd.DataFrame(columns=['index', 'angle', 'height', 'prominence'])
        
        peaks_df = pd.DataFrame({
            'index': peak_indices,
            'angle': angles[peak_indices],
            'height': signal[peak_indices],
            'prominence': properties['prominences']
        })
    else:
        # Fallback: simple local maxima detection
        peak_indices = []
        prominences = []
        
        for i in range(1, len(signal) - 1):
            if signal[i] > signal[i-1] and signal[i] > signal[i+1] and signal[i] > 0:
                # Estimate prominence as height above local minimum
                local_window = 20
                start = max(0, i - local_window)
                end = min(len(signal), i + local_window + 1)
                local_min = np.min(signal[start:end])
                prominence = signal[i] - local_min
                
                if prominence >= prominence_threshold:
                    peak_indices.append(i)
                    prominences.append(prominence)
        
        if len(peak_indices) == 0:
            return pd.DataFrame(columns=['index', 'angle', 'height', 'prominence'])
        
        peaks_df = pd.DataFrame({
            'index': peak_indices,
            'angle': angles[peak_indices],
            'height': signal[peak_indices],
            'prominence': prominences
        })
    
    # Sort by prominence (descending)
    peaks_df = peaks_df.sort_values('prominence', ascending=False).reset_index(drop=True)
    
    return peaks_df


def choose_small_peak(
    peaks: pd.DataFrame,
    min_separation_deg: float = Config.MIN_SEPARATION_DEG,
    min_height_ratio: float = Config.MIN_HEIGHT_RATIO
) -> Tuple[Optional[Dict], Optional[Dict], str]:
    """
    Identify the main peak and the small direct-light peak.
    
    Args:
        peaks: DataFrame of detected peaks
        min_separation_deg: Minimum angle separation between main and small peak
        min_height_ratio: Small peak must be at least this fraction of main peak height
        
    Returns:
        Tuple of (main_peak_dict, small_peak_dict, notes_string)
        Either peak dict can be None if not found
    """
    if len(peaks) == 0:
        return None, None, "No peaks detected"
    
    # Main peak is the one with largest prominence
    main_peak = peaks.iloc[0].to_dict()
    
    if len(peaks) == 1:
        return main_peak, None, "Only one peak detected"
    
    # Look for small peak: must be separated from main peak and above noise
    min_height = main_peak['height'] * min_height_ratio
    
    candidates = []
    for idx, peak in peaks.iloc[1:].iterrows():
        angle_separation = abs(peak['angle'] - main_peak['angle'])
        
        if angle_separation >= min_separation_deg and peak['height'] >= min_height:
            candidates.append(peak)
    
    if len(candidates) == 0:
        return main_peak, None, f"No secondary peak found (separation > {min_separation_deg}°, height > {min_height_ratio:.1%} of main)"
    
    # Choose the most prominent candidate
    small_peak = candidates[0].to_dict()
    notes = f"Small peak found at {small_peak['angle']:.2f}° (prominence: {small_peak['prominence']:.4f})"
    
    return main_peak, small_peak, notes


# ============================================================================
# PLOTTING
# ============================================================================

def plot_scan(
    file_name: str,
    angles: np.ndarray,
    raw_intensity: np.ndarray,
    smoothed: np.ndarray,
    baseline_corrected: np.ndarray,
    main_peak: Optional[Dict],
    small_peak: Optional[Dict],
    output_path: Path,
    debug: bool = False
) -> None:
    """
    Create a diagnostic plot for a scan.
    
    Args:
        file_name: Name of the data file
        angles: Angle values
        raw_intensity: Raw intensity data
        smoothed: Smoothed intensity
        baseline_corrected: Baseline-corrected intensity
        main_peak: Main peak dictionary (or None)
        small_peak: Small peak dictionary (or None)
        output_path: Path to save the plot
        debug: If True, add extra information to plot
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    # Top panel: Raw and smoothed
    ax1.plot(angles, raw_intensity, 'k-', alpha=0.3, linewidth=0.5, label='Raw')
    ax1.plot(angles, smoothed, 'b-', linewidth=1.5, label='Smoothed')
    ax1.set_ylabel('Intensity (V)')
    ax1.set_title(f'Step 1 Analysis: {file_name}')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)
    
    # Bottom panel: Baseline-corrected with peaks
    ax2.plot(angles, baseline_corrected, 'g-', linewidth=2, label='Baseline-corrected')
    ax2.axhline(0, color='k', linestyle='--', alpha=0.3)
    
    # Mark peaks
    if main_peak:
        ax2.axvline(main_peak['angle'], color='red', linestyle='--', linewidth=2, alpha=0.7, label='Main peak')
        ax2.plot(main_peak['angle'], main_peak['height'], 'ro', markersize=10)
    
    if small_peak:
        ax2.axvline(small_peak['angle'], color='orange', linestyle='--', linewidth=2, alpha=0.7, label='Small peak (θ_init)')
        ax2.plot(small_peak['angle'], small_peak['height'], 'o', color='orange', markersize=10)
        
        # Annotate small peak angle
        ax2.annotate(
            f"θ_init = {small_peak['angle']:.2f}°",
            xy=(small_peak['angle'], small_peak['height']),
            xytext=(10, 20),
            textcoords='offset points',
            fontsize=12,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', lw=2)
        )
    
    ax2.set_xlabel('Angle (degrees)')
    ax2.set_ylabel('Intensity (V, baseline-corrected)')
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================

def analyze_step1_file(
    file_path: Path,
    config: Config,
    output_dir: Path
) -> Dict:
    """
    Analyze a single Step 1 file.
    
    Args:
        file_path: Path to the data file
        config: Configuration object
        output_dir: Directory for output files
        
    Returns:
        Dictionary with analysis results
    """
    result = {
        'file': file_path.name,
        'theta_init_deg': np.nan,
        'theta_init_method': 'auto_peak_detection',
        'notes': ''
    }
    
    try:
        # Load data
        df = load_scan(file_path)
        
        # Preprocess
        angles, raw, smoothed, baseline_corrected = preprocess_scan(
            df,
            smooth_window=config.SMOOTH_WINDOW,
            smooth_poly=config.SMOOTH_POLY,
            baseline_percentile=config.BASELINE_PERCENTILE
        )
        
        # Detect peaks
        peaks = detect_peaks(
            angles,
            baseline_corrected,
            prominence_threshold=config.PROMINENCE_THRESHOLD
        )
        
        # Choose main and small peaks
        main_peak, small_peak, notes = choose_small_peak(
            peaks,
            min_separation_deg=config.MIN_SEPARATION_DEG,
            min_height_ratio=config.MIN_HEIGHT_RATIO
        )
        
        result['notes'] = notes
        
        if small_peak:
            result['theta_init_deg'] = small_peak['angle']
        else:
            warnings.warn(f"Could not find small peak in {file_path.name}: {notes}")
        
        # Create plot
        plot_path = output_dir / 'plots' / f"{file_path.stem}_step1.png"
        plot_scan(
            file_path.name,
            angles,
            raw,
            smoothed,
            baseline_corrected,
            main_peak,
            small_peak,
            plot_path,
            debug=config.DEBUG
        )
        
        if config.DEBUG and len(peaks) > 0:
            # Save peak table for debugging
            peaks_path = output_dir / 'debug' / f"{file_path.stem}_peaks.csv"
            peaks_path.parent.mkdir(parents=True, exist_ok=True)
            peaks.to_csv(peaks_path, index=False)
        
    except Exception as e:
        result['notes'] = f"ERROR: {e}"
        warnings.warn(f"Failed to analyze {file_path.name}: {e}")
    
    return result


def compute_statistics(theta_init_values: List[float]) -> Dict[str, float]:
    """
    Compute statistical measures for theta_init values.
    
    Args:
        theta_init_values: List of theta_init measurements
        
    Returns:
        Dictionary with statistical measures
    """
    values = np.array([v for v in theta_init_values if not np.isnan(v)])
    
    if len(values) == 0:
        return {
            'n': 0,
            'mean': np.nan,
            'std': np.nan,
            'sem': np.nan,
            'min': np.nan,
            'max': np.nan,
            'half_range': np.nan
        }
    
    stats = {
        'n': len(values),
        'mean': np.mean(values),
        'std': np.std(values, ddof=1) if len(values) > 1 else 0.0,
        'min': np.min(values),
        'max': np.max(values)
    }
    
    stats['sem'] = stats['std'] / np.sqrt(stats['n']) if stats['n'] > 0 else np.nan
    stats['half_range'] = (stats['max'] - stats['min']) / 2.0
    
    return stats


def main():
    """Main entry point for Step 1 analysis."""
    parser = argparse.ArgumentParser(
        description='Analyze Blackbody Lab Step 1 data to extract theta_init',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input/output
    parser.add_argument('--data_dir', type=str, default=Config.DATA_DIR,
                        help='Directory containing data files')
    parser.add_argument('--pattern', type=str, default=Config.PATTERN,
                        help='Glob pattern for Step 1 files')
    parser.add_argument('--out_dir', type=str, default=Config.OUTPUT_DIR,
                        help='Output directory')
    
    # Processing parameters
    parser.add_argument('--prominence', type=float, default=Config.PROMINENCE_THRESHOLD,
                        help='Minimum peak prominence (V)')
    parser.add_argument('--min_separation_deg', type=float, default=Config.MIN_SEPARATION_DEG,
                        help='Minimum angle separation between peaks (degrees)')
    parser.add_argument('--smooth_window', type=int, default=Config.SMOOTH_WINDOW,
                        help='Smoothing window size')
    parser.add_argument('--smooth_poly', type=int, default=Config.SMOOTH_POLY,
                        help='Polynomial order for Savitzky-Golay')
    
    # Debug
    parser.add_argument('--debug', action='store_true',
                        help='Enable debug mode (save extra files)')
    
    args = parser.parse_args()
    
    # Update config from args
    config = Config()
    config.DATA_DIR = args.data_dir
    config.PATTERN = args.pattern
    config.OUTPUT_DIR = args.out_dir
    config.PROMINENCE_THRESHOLD = args.prominence
    config.MIN_SEPARATION_DEG = args.min_separation_deg
    config.SMOOTH_WINDOW = args.smooth_window
    config.SMOOTH_POLY = args.smooth_poly
    config.DEBUG = args.debug
    
    # Find Step 1 files
    data_dir = Path(config.DATA_DIR)
    if not data_dir.exists():
        print(f"ERROR: Data directory not found: {data_dir}")
        return
    
    step1_files = sorted(data_dir.glob(config.PATTERN))
    
    if len(step1_files) == 0:
        print(f"ERROR: No Step 1 files found matching pattern: {config.PATTERN}")
        return
    
    print(f"Found {len(step1_files)} Step 1 files")
    print("=" * 70)
    
    # Create output directory
    output_dir = Path(config.OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Analyze each file
    results = []
    for file_path in step1_files:
        print(f"Processing: {file_path.name}")
        result = analyze_step1_file(file_path, config, output_dir)
        results.append(result)
        
        if not np.isnan(result['theta_init_deg']):
            print(f"  → θ_init = {result['theta_init_deg']:.3f}°")
        else:
            print(f"  → Could not determine θ_init: {result['notes']}")
    
    # Save results table
    results_df = pd.DataFrame(results)
    table_path = output_dir / 'step1_theta_init_table.csv'
    results_df.to_csv(table_path, index=False)
    print(f"\nResults table saved to: {table_path}")
    
    # Compute statistics
    theta_init_values = results_df['theta_init_deg'].values
    stats = compute_statistics(theta_init_values)
    
    # Generate summary
    summary_lines = []
    summary_lines.append("=" * 70)
    summary_lines.append("STEP 1 ANALYSIS SUMMARY: theta_init Extraction")
    summary_lines.append("=" * 70)
    summary_lines.append(f"Files processed: {len(results)}")
    summary_lines.append(f"Successful extractions: {stats['n']}")
    summary_lines.append("")
    
    if stats['n'] > 0:
        summary_lines.append(f"Mean θ_init:           {stats['mean']:.4f}°")
        summary_lines.append(f"Sample std dev (σ):    {stats['std']:.4f}°")
        summary_lines.append(f"Std error of mean:     {stats['sem']:.4f}°")
        summary_lines.append(f"Range:                 [{stats['min']:.4f}°, {stats['max']:.4f}°]")
        summary_lines.append(f"Half-range:            {stats['half_range']:.4f}°")
        summary_lines.append("")
        summary_lines.append("UNCERTAINTY RECOMMENDATIONS:")
        summary_lines.append(f"  - For measurement scatter:     ± {stats['std']:.4f}° (sample std dev)")
        summary_lines.append(f"  - For uncertainty in mean:     ± {stats['sem']:.4f}° (std error)")
        summary_lines.append(f"  - Alternative (half-range):    ± {stats['half_range']:.4f}°")
        summary_lines.append("")
        summary_lines.append("Recommended: Use sample std dev (σ) for typical measurement uncertainty,")
        summary_lines.append("or std error of mean (SEM) if reporting the mean value specifically.")
    else:
        summary_lines.append("No successful theta_init extractions!")
    
    summary_lines.append("=" * 70)
    
    summary_text = "\n".join(summary_lines)
    
    # Save summary
    summary_path = output_dir / 'step1_summary.txt'
    with open(summary_path, 'w') as f:
        f.write(summary_text)
    
    # Print summary to terminal
    print("\n" + summary_text)
    print(f"\nSummary saved to: {summary_path}")
    print(f"Plots saved to: {output_dir / 'plots'}/")
    
    print("\n✓ Analysis complete!")


if __name__ == '__main__':
    main()
