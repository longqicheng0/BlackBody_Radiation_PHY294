# Step 1 Analysis - Quick Reference

## What This Does
Automatically finds the small "direct light" peak (θ_init) in each Step 1 blackbody scan.

## Running the Analysis

### Option 1: Quick Start (Recommended)
```bash
./run_step1.sh
```

### Option 2: Manual Run
```bash
# Install dependencies first (one time only)
pip3 install -r requirements.txt

# Run analysis
python3 -m bb_step1
```

## Output Files

After running, check the `outputs/` directory:

- **`outputs/plots/`** - Individual scan plots with peak annotations
- **`outputs/step1_theta_init_table.csv`** - Table of extracted θ_init values
- **`outputs/step1_summary.txt`** - Statistical summary with uncertainties

## Customization

### Adjust peak detection sensitivity:
```bash
python3 -m bb_step1 --prominence 0.005  # More sensitive (finds smaller peaks)
python3 -m bb_step1 --prominence 0.02   # Less sensitive (only large peaks)
```

### Change angle separation threshold:
```bash
python3 -m bb_step1 --min_separation_deg 3.0  # Require larger separation
```

### Enable debug mode:
```bash
python3 -m bb_step1 --debug  # Saves extra diagnostic files
```

### Use different data directory:
```bash
python3 -m bb_step1 --data_dir path/to/data --pattern "step1*.txt"
```

## Interpreting Results

The summary includes three uncertainty measures:

1. **Sample std dev (σ)** - Scatter in measurements (recommended for typical uncertainty)
2. **Std error (SEM)** - Uncertainty in the mean value (use if reporting mean)
3. **Half-range** - Conservative estimate: (max - min) / 2

## Troubleshooting

**Problem:** No peaks detected
- **Solution:** Lower prominence threshold: `--prominence 0.005`

**Problem:** Wrong peak selected
- **Solution:** Increase separation: `--min_separation_deg 3.0`
- **Solution:** Check plots to understand peak structure

**Problem:** Import errors
- **Solution:** Reinstall dependencies: `pip3 install -r requirements.txt`

## Peak Selection Algorithm

The script:
1. Detects ALL peaks in the scan
2. Identifies the MAIN peak (largest prominence = blackbody emission)
3. Finds the SMALL peak (θ_init) that is:
   - Separated from main peak by > 2° (configurable)
   - Height > 5% of main peak
   - Most prominent among candidates

## Data Format Requirements

Step 1 files should be:
- Tab or space-separated text
- Two columns: angle (degrees), intensity (volts)
- Headers automatically skipped
- Named matching pattern: `step1*.txt`

## Questions?

See the full README.md for detailed documentation.
