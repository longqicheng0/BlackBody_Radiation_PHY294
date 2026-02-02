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
from typing import Dict, List, Optional, Tuple, Sequence

import matplotlib
matplotlib.use('Agg')  # Non-GUI backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Optional scipy for peak detection
try:
    from scipy.signal import find_peaks
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
    PATTERN = "**/[sS]tep1*.txt"
    OUTPUT_DIR = "outputs"
    
    # Preprocessing
    # (Baseline correction removed; use raw data directly.)
    
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
    df: pd.DataFrame
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Preprocess scan data (no baseline correction).
    
    Args:
        df: DataFrame with 'angle' and 'intensity' columns
        
    Returns:
        Tuple of (angles, raw_intensity)
    """
    angles = df['angle'].values
    raw_intensity = df['intensity'].values
    
    return angles, raw_intensity


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
        signal: Intensity signal
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
            height=0  # Only positive peaks
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


def choose_theta_init_peak(
    peaks: pd.DataFrame,
    min_separation_deg: float = Config.MIN_SEPARATION_DEG,
    min_height_ratio: float = Config.MIN_HEIGHT_RATIO
) -> Tuple[Optional[Dict], Optional[Dict], str]:
    """
    Identify the theta_init peak (highest intensity) and the secondary peak.
    
    Args:
        peaks: DataFrame of detected peaks
        min_separation_deg: Minimum angle separation between peaks
        min_height_ratio: Secondary peak must be at least this fraction of theta_init height
        
    Returns:
        Tuple of (theta_init_peak_dict, secondary_peak_dict, notes_string)
        Either peak dict can be None if not found
    """
    if len(peaks) == 0:
        return None, None, "No peaks detected"
    
    # Theta_init is the peak with the highest prominence (largest intensity feature)
    theta_init_peak = peaks.iloc[0].to_dict()
    
    if len(peaks) == 1:
        return theta_init_peak, None, "Only one peak detected"
    
    # Look for secondary peak: separated from theta_init and above noise
    min_height = theta_init_peak['height'] * min_height_ratio
    
    candidates = []
    for idx, peak in peaks.iloc[1:].iterrows():
        angle_separation = abs(peak['angle'] - theta_init_peak['angle'])
        
        if angle_separation >= min_separation_deg and peak['height'] >= min_height:
            candidates.append(peak)
    
    if len(candidates) == 0:
        return theta_init_peak, None, f"No secondary peak found (separation > {min_separation_deg}°, height > {min_height_ratio:.1%} of theta_init)"
    
    # Choose the most prominent candidate
    secondary_peak = candidates[0].to_dict()
    notes = f"Secondary peak found at {secondary_peak['angle']:.2f}° (prominence: {secondary_peak['prominence']:.4f})"
    
    return theta_init_peak, secondary_peak, notes


# ============================================================================
# PLOTTING
# ============================================================================

def plot_scan(
    file_name: str,
    angles: np.ndarray,
    raw_intensity: np.ndarray,
    theta_init_peak: Optional[Dict],
    secondary_peak: Optional[Dict],
    output_path: Path,
    debug: bool = False,
    uncertainty_deg: Optional[float] = None
) -> None:
    """
    Create a diagnostic plot for a scan.
    
    Args:
        file_name: Name of the data file
        angles: Angle values
        raw_intensity: Raw intensity data
        theta_init_peak: Theta_init peak dictionary (or None)
        secondary_peak: Secondary peak dictionary (or None)
        output_path: Path to save the plot
        debug: If True, add extra information to plot
        uncertainty_deg: Uncertainty in theta_init (optional, for annotation)
    """
    fig, ax2 = plt.subplots(1, 1, figsize=(12, 6))
    
    # Single panel: Raw data with peaks
    ax2.plot(angles, raw_intensity, 'g-', linewidth=2, label='Raw')
    ax2.set_title(f'Step 1 Analysis: {file_name}')
    
    # Mark peaks
    if theta_init_peak:
        ax2.axvline(theta_init_peak['angle'], color='orange', linestyle='--', linewidth=2, alpha=0.7, label='θ_init peak')
        ax2.plot(theta_init_peak['angle'], theta_init_peak['height'], 'o', color='orange', markersize=10)
    
    if secondary_peak:
        ax2.axvline(secondary_peak['angle'], color='red', linestyle='--', linewidth=2, alpha=0.7, label='Secondary peak')
        ax2.plot(secondary_peak['angle'], secondary_peak['height'], 'ro', markersize=10)
        
    # Annotate theta_init angle (with uncertainty if provided)
    if theta_init_peak:
        if uncertainty_deg is not None:
            annotation_text = f"θ_init = {theta_init_peak['angle']:.2f}° ± {uncertainty_deg:.4f}°"
        else:
            annotation_text = f"θ_init = {theta_init_peak['angle']:.2f}°"
        
        ax2.annotate(
            annotation_text,
            xy=(theta_init_peak['angle'], theta_init_peak['height']),
            xytext=(10, 20),
            textcoords='offset points',
            fontsize=12,
            fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', lw=2)
        )
    
    ax2.set_xlabel('Angle (degrees)')
    ax2.set_ylabel('Intensity (V)')
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
        'theta_init_method': 'auto_peak_detection_highest_peak',
        'notes': ''
    }
    
    try:
        # Load data
        df = load_scan(file_path)
        
        # Preprocess
        angles, raw = preprocess_scan(df)
        
        # Detect peaks
        peaks = detect_peaks(
            angles,
            raw,
            prominence_threshold=config.PROMINENCE_THRESHOLD
        )
        
        # Choose theta_init and secondary peaks
        theta_init_peak, secondary_peak, notes = choose_theta_init_peak(
            peaks,
            min_separation_deg=config.MIN_SEPARATION_DEG,
            min_height_ratio=config.MIN_HEIGHT_RATIO
        )
        
        result['notes'] = notes
        
        if theta_init_peak:
            result['theta_init_deg'] = theta_init_peak['angle']
        else:
            warnings.warn(f"Could not find small peak in {file_path.name}: {notes}")
        
        # Create plot
        plot_path = output_dir / 'plots' / f"{file_path.stem}_step1.png"
        plot_scan(
            file_path.name,
            angles,
            raw,
            theta_init_peak,
            secondary_peak,
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


def compute_step1_theta_init_uncertainty(theta_init_deg: Sequence[float]) -> dict:
    """
    Compute uncertainty statistics for Step 1 theta_init across repeated measurements.
    
    Uses standard 2nd-year physics definitions:
      - sample standard deviation (N-1 in denominator)
      - standard error of the mean (SEM = sd / sqrt(N))
    
    Args:
        theta_init_deg: Sequence of theta_init measurements in degrees
        
    Returns:
        Dictionary with keys:
          'n', 'mean_deg', 'sd_deg', 'sem_deg', 'min_deg', 'max_deg', 'half_range_deg'
          
    Raises:
        ValueError: If n < 2 or any NaN exists
    """
    # Check for NaNs
    values = np.array(theta_init_deg)
    nan_mask = np.isnan(values)
    if nan_mask.any():
        nan_indices = np.where(nan_mask)[0]
        raise ValueError(f"Found NaN theta_init values at indices: {nan_indices.tolist()}")
    
    n = len(values)
    if n < 2:
        raise ValueError(f"Need at least 2 measurements for uncertainty analysis, got {n}")
    
    # Standard 2nd-year physics error analysis:
    # Sample standard deviation (N-1), SEM = sd/sqrt(N)
    mean_deg = float(np.mean(values))
    sd_deg = float(np.std(values, ddof=1))  # ddof=1 for sample std dev
    sem_deg = sd_deg / np.sqrt(n)
    min_deg = float(np.min(values))
    max_deg = float(np.max(values))
    half_range_deg = (max_deg - min_deg) / 2.0
    
    return {
        'n': n,
        'mean_deg': mean_deg,
        'sd_deg': sd_deg,
        'sem_deg': sem_deg,
        'min_deg': min_deg,
        'max_deg': max_deg,
        'half_range_deg': half_range_deg
    }


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


def plot_summary(results_df: pd.DataFrame, output_dir: Path) -> None:
    """
    Create summary plots and tables for theta_init results.

    Args:
        results_df: DataFrame with theta_init results
        output_dir: Base output directory
    """
    summary_dir = output_dir / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    # Clean values
    clean_df = results_df.dropna(subset=['theta_init_deg']).copy()
    if clean_df.empty:
        return

    # Summary table
    values = clean_df['theta_init_deg'].values
    summary_stats = {
        'count': len(values),
        'mean': float(np.mean(values)),
        'std': float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
        'sem': float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else 0.0,
        'min': float(np.min(values)),
        'max': float(np.max(values)),
        'median': float(np.median(values)),
        'q1': float(np.percentile(values, 25)),
        'q3': float(np.percentile(values, 75)),
        'iqr': float(np.percentile(values, 75) - np.percentile(values, 25))
    }
    summary_table = pd.DataFrame([summary_stats])
    summary_table.to_csv(summary_dir / 'step1_theta_init_summary.csv', index=False)

    # Plot 1: theta_init by file
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(clean_df))
    ax.plot(x, values, 'o-', color='tab:blue', label='θ_init')
    ax.axhline(summary_stats['mean'], color='tab:orange', linestyle='--', label='Mean')
    ax.set_xticks(x)
    ax.set_xticklabels(clean_df['file'].tolist(), rotation=45, ha='right')
    ax.set_ylabel('θ_init (degrees)')
    ax.set_title('θ_init by File')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig(summary_dir / 'step1_theta_init_by_file.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Plot 2: histogram
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(values, bins='auto', color='tab:green', alpha=0.7, edgecolor='black')
    ax.axvline(summary_stats['mean'], color='tab:orange', linestyle='--', label='Mean')
    ax.set_xlabel('θ_init (degrees)')
    ax.set_ylabel('Count')
    ax.set_title('θ_init Distribution')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig(summary_dir / 'step1_theta_init_hist.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Plot 3: boxplot
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.boxplot(values, vert=True, patch_artist=True,
               boxprops=dict(facecolor='lightblue', color='black'),
               medianprops=dict(color='red'))
    ax.set_ylabel('θ_init (degrees)')
    ax.set_title('θ_init Boxplot')
    ax.grid(True, axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(summary_dir / 'step1_theta_init_box.png', dpi=150, bbox_inches='tight')
    plt.close()


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
    
    # Compute uncertainty statistics and add to each result row
    results_df = pd.DataFrame(results)
    theta_init_values = results_df['theta_init_deg'].values
    
    # Check if we have the expected number of files
    if len(step1_files) != 5:
        warnings.warn(f"Expected 5 Step 1 files, found {len(step1_files)}")
    
    # Compute uncertainties
    try:
        uncertainty_stats = compute_step1_theta_init_uncertainty(theta_init_values)
        
        # Add uncertainty columns to each row (same values repeated)
        results_df['theta_init_mean_deg'] = uncertainty_stats['mean_deg']
        results_df['theta_init_sd_deg'] = uncertainty_stats['sd_deg']
        results_df['theta_init_sem_deg'] = uncertainty_stats['sem_deg']
        results_df['theta_init_n'] = uncertainty_stats['n']
        
        # Print uncertainty summary to terminal
        print(f"\n{'='*70}")
        print("UNCERTAINTY ANALYSIS (Step 1 θ_init)")
        print(f"{'='*70}")
        print(f"N (measurements):          {uncertainty_stats['n']}")
        print(f"Mean:                      {uncertainty_stats['mean_deg']:.4f}°")
        print(f"Sample std dev (σ):        {uncertainty_stats['sd_deg']:.4f}°")
        print(f"Std error of mean (SEM):   {uncertainty_stats['sem_deg']:.4f}°")
        print(f"Range:                     [{uncertainty_stats['min_deg']:.4f}°, {uncertainty_stats['max_deg']:.4f}°]")
        print(f"Half-range:                {uncertainty_stats['half_range_deg']:.4f}°")
        print(f"{'='*70}\n")
        
        # Regenerate plots with uncertainty annotations
        print("Updating plots with uncertainty annotations...")
        uncertainty_deg = uncertainty_stats['sd_deg']  # Use sample std dev for uncertainty
        for i, (file_path, result) in enumerate(zip(step1_files, results)):
            if not np.isnan(result['theta_init_deg']):
                # Re-analyze to get peak data for plotting
                try:
                    df = load_scan(file_path)
                    angles, raw = preprocess_scan(df)
                    peaks = detect_peaks(angles, raw, prominence_threshold=config.PROMINENCE_THRESHOLD)
                    theta_init_peak, secondary_peak, _ = choose_theta_init_peak(
                        peaks, min_separation_deg=config.MIN_SEPARATION_DEG, min_height_ratio=config.MIN_HEIGHT_RATIO
                    )
                    plot_path = output_dir / 'plots' / f"{file_path.stem}_step1.png"
                    plot_scan(file_path.name, angles, raw, theta_init_peak, secondary_peak, plot_path, 
                             debug=config.DEBUG, uncertainty_deg=uncertainty_deg)
                except Exception as e:
                    warnings.warn(f"Could not update plot for {file_path.name}: {e}")
        
    except ValueError as e:
        print(f"\nWARNING: Could not compute uncertainties: {e}")
        # Add NaN columns if computation failed
        results_df['theta_init_mean_deg'] = np.nan
        results_df['theta_init_sd_deg'] = np.nan
        results_df['theta_init_sem_deg'] = np.nan
        results_df['theta_init_n'] = 0
    
    # Save results table with uncertainty columns
    table_path = output_dir / 'step1_theta_init_table.csv'
    results_df.to_csv(table_path, index=False)
    print(f"Results table saved to: {table_path}")

    # Summary plots and table
    plot_summary(results_df, output_dir)
    
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
    print(f"Summary plots saved to: {output_dir / 'summary'}/")
    
    print("\n✓ Analysis complete!")


if __name__ == '__main__':
    main()
