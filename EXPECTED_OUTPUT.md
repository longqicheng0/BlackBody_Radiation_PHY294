# Expected Output Structure

After running the Step 1 analysis, you should see this structure:

```
outputs/
├── plots/
│   ├── Step1_1_step1.png
│   ├── step1_2_step1.png
│   ├── step1_3(def)_step1.png
│   ├── step1_4_step1.png
│   ├── step1_5_step1.png
│   └── step1_6_step1.png
├── step1_theta_init_table.csv
└── step1_summary.txt
```

## File Descriptions

### Individual Plots (`plots/*.png`)
Each plot contains two panels:

**Top panel:**
- Raw intensity data (thin black line)
- Smoothed intensity (thick blue line)

**Bottom panel:**
- Baseline-corrected signal (green line)
- Main blackbody peak (red dashed line)
- Small direct-light peak / θ_init (orange dashed line)
- Annotation showing θ_init value

### Results Table (`step1_theta_init_table.csv`)
CSV format with columns:
```csv
file,theta_init_deg,theta_init_method,notes
Step1_1.txt,42.35,auto_peak_detection,"Small peak found at 42.35° (prominence: 0.0234)"
step1_2.txt,42.41,auto_peak_detection,"Small peak found at 42.41° (prominence: 0.0198)"
...
```

### Statistical Summary (`step1_summary.txt`)
Text file containing:
```
======================================================================
STEP 1 ANALYSIS SUMMARY: theta_init Extraction
======================================================================
Files processed: 6
Successful extractions: 6

Mean θ_init:           42.3783°
Sample std dev (σ):    0.0523°
Std error of mean:     0.0213°
Range:                 [42.3100°, 42.4500°]
Half-range:            0.0700°

UNCERTAINTY RECOMMENDATIONS:
  - For measurement scatter:     ± 0.0523° (sample std dev)
  - For uncertainty in mean:     ± 0.0213° (std error)
  - Alternative (half-range):    ± 0.0700°

Recommended: Use sample std dev (σ) for typical measurement uncertainty,
or std error of mean (SEM) if reporting the mean value specifically.
======================================================================
```

## Sample Terminal Output

When you run `./run_step1.sh` or `python3 -m bb_step1`, expect:

```
Found 6 Step 1 files
======================================================================
Processing: Step1_1.txt
  → θ_init = 42.350°
Processing: step1_2.txt
  → θ_init = 42.410°
Processing: step1_3(def).txt
  → θ_init = 42.385°
Processing: step1_4.txt
  → θ_init = 42.340°
Processing: step1_5.txt
  → θ_init = 42.450°
Processing: step1_6.txt
  → θ_init = 42.310°

Results table saved to: outputs/step1_theta_init_table.csv

======================================================================
STEP 1 ANALYSIS SUMMARY: theta_init Extraction
======================================================================
Files processed: 6
Successful extractions: 6

Mean θ_init:           42.3783°
Sample std dev (σ):    0.0523°
Std error of mean:     0.0213°
Range:                 [42.3100°, 42.4500°]
Half-range:            0.0700°

UNCERTAINTY RECOMMENDATIONS:
  - For measurement scatter:     ± 0.0523° (sample std dev)
  - For uncertainty in mean:     ± 0.0213° (std error)
  - Alternative (half-range):    ± 0.0700°

Recommended: Use sample std dev (σ) for typical measurement uncertainty,
or std error of mean (SEM) if reporting the mean value specifically.
======================================================================

Summary saved to: outputs/step1_summary.txt
Plots saved to: outputs/plots/

✓ Analysis complete!
```

## Debug Mode Output

When running with `--debug` flag:

```
outputs/
├── debug/
│   ├── Step1_1_peaks.csv      # All detected peaks with metadata
│   ├── step1_2_peaks.csv
│   └── ...
├── plots/
│   └── ...
├── step1_theta_init_table.csv
└── step1_summary.txt
```

The debug peak tables show all detected peaks:
```csv
index,angle,height,prominence
125,42.35,0.045,0.0234
278,65.80,0.234,0.1856
```

## Troubleshooting

If you see **warnings** like:
```
WARNING: Could not find small peak in step1_X.txt: No secondary peak found...
```

This means the algorithm couldn't confidently identify θ_init. Check:
1. The corresponding plot to see peak structure
2. Try adjusting `--prominence` or `--min_separation_deg`
3. Verify the data file has expected format

If **no outputs** are generated:
- Check that data files exist in `Blackbody_Lab_Data/`
- Ensure dependencies are installed: `pip3 install -r requirements.txt`
- Try running with `--debug` for more information
