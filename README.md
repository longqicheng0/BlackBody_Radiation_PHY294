# BlackBody Radiation Lab - PHY294

Python analysis tools for processing blackbody radiation lab data.

## Overview

This repository contains automated analysis tools for PHY294 blackbody radiation experiments. Currently implements **Step 1 analysis** to extract the direct/unrefracted light peak angle (θ_init).

## Repository Structure

```
.
├── Blackbody_Lab_Data/          # Raw experimental data files (.txt)
│   ├── step1_*.txt              # Step 1 scans (angle vs intensity)
│   └── step2_*.txt              # Step 2 scans (voltage-dependent)
├── bb_step1/                    # Step 1 analysis package
│   ├── __init__.py
│   └── __main__.py              # Main analysis script
├── outputs/                     # Generated analysis outputs
│   ├── plots/                   # Individual scan plots
│   ├── step1_theta_init_table.csv
│   └── step1_summary.txt
├── requirements.txt             # Python dependencies
├── README.md                    # This file
└── LICENSE
```

## Installation

### Prerequisites
- Python 3.9 or higher (3.10+ recommended)
- pip

### Setup

1. Clone or navigate to this repository:
   ```bash
   cd BlackBody_Radiation_PHY294
   ```

2. Install dependencies:
   ```bash
   pip3 install -r requirements.txt
   ```
   
   Or if you prefer using a virtual environment:
   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On macOS/Linux
   pip install -r requirements.txt
   ```

## Usage

### Quick Start

The easiest way to run the analysis:

```bash
./run_step1.sh
```

This script will automatically:
- Install required dependencies
- Run the Step 1 analysis on all data files
- Generate plots and statistical summaries

### Step 1 Analysis: θ_init Extraction

The Step 1 analysis identifies the small "direct light" peak that occurs when light bypasses the prism, extracting the angle θ_init for each scan.

**Run the analysis manually:**
```bash
python3 -m bb_step1
```

**Custom options:**
```bash
python3 -m bb_step1 \
  --data_dir Blackbody_Lab_Data \
  --pattern "**/step1*.txt" \
  --out_dir outputs \
  --prominence 0.01 \
  --min_separation_deg 2.0 \
  --smooth_window 11 \
  --debug
```

**Key parameters:**
- `--data_dir`: Directory containing data files (default: `Blackbody_Lab_Data`)
- `--pattern`: Glob pattern for Step 1 files (default: `**/step1*.txt`)
- `--out_dir`: Output directory (default: `outputs`)
- `--prominence`: Minimum peak prominence in volts (default: 0.01)
- `--min_separation_deg`: Minimum angle separation between main and small peaks (default: 2.0°)
- `--smooth_window`: Smoothing window size (default: 11)
- `--debug`: Enable debug mode (saves extra diagnostic files)

### Outputs

After running the analysis, the following files are generated:

1. **`outputs/plots/<filename>_step1.png`** - Individual scan plots showing:
   - Raw and smoothed intensity curves
   - Baseline-corrected signal
   - Main blackbody peak (red line)
   - Small direct-light peak θ_init (orange line, annotated)

2. **`outputs/step1_theta_init_table.csv`** - Table with columns:
   - `file`: Data file name
   - `theta_init_deg`: Extracted θ_init angle
   - `theta_init_method`: Method used (auto_peak_detection)
   - `notes`: Additional information or error messages

3. **`outputs/step1_summary.txt`** - Statistical summary including:
   - Mean θ_init
   - Sample standard deviation (measurement scatter)
   - Standard error of the mean
   - Range and half-range
   - Uncertainty recommendations

### How It Works

The analysis pipeline:

1. **Data Loading**: Automatically discovers and loads Step 1 files (tab-separated, skips headers)

2. **Preprocessing**:
   - Sorts data by angle
   - Applies smoothing (Savitzky-Golay filter if scipy available, else moving average)
   - Baseline correction (subtracts low-percentile baseline)

3. **Peak Detection**:
   - Identifies all peaks above a prominence threshold
   - Classifies the **main peak** (largest prominence - the blackbody emission)
   - Identifies the **small peak** (θ_init) as the most prominent secondary peak that:
     - Is separated from the main peak by > 2° (configurable)
     - Has sufficient height (> 5% of main peak)
     - Is not noise

4. **Output Generation**:
   - Creates diagnostic plots with peak annotations
   - Compiles results table
   - Computes statistics and uncertainty estimates

### Uncertainty Recommendations

The analysis reports three uncertainty measures:

- **Sample standard deviation (σ)**: Recommended for typical measurement uncertainty/scatter
- **Standard error of mean (SEM = σ/√n)**: Use when reporting uncertainty in the mean value
- **Half-range**: Alternative conservative estimate

## Data Format

Step 1 data files (`step1_*.txt`) should be tab or space-separated text with:
- Optional header lines (automatically skipped)
- Two numeric columns: angle (degrees) and intensity (volts)
- Example:
  ```
  Sensor Position (Degrees)    Light Intensity (Volts)
  10.5    0.023
  11.0    0.025
  ...
  ```

## Troubleshooting

**No peaks detected:**
- Try reducing `--prominence` (e.g., `--prominence 0.005`)
- Check plots to verify data quality
- Enable `--debug` mode to inspect peak detection tables

**Wrong peak selected:**
- Adjust `--min_separation_deg` to enforce better separation
- Check if baseline correction is working (view plots)

**Import errors:**
- Ensure all dependencies are installed: `pip install -r requirements.txt`
- scipy is optional; fallback methods are available

## License

See [LICENSE](LICENSE) file for details.
