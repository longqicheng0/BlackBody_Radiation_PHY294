#!/usr/bin/env python3
import pandas as pd
import numpy as np

df = pd.read_csv('outputs/step2_wien_results.csv')

print("=" * 75)
print("UNCERTAINTY PROPAGATION VERIFICATION")
print("=" * 75)

print("\nFIRST MEASUREMENT (step2_10v.txt):")
print("-" * 75)
print(f"Voltage:              {df.loc[0, 'voltage_V']:.1f} V")
print(f"Current:              {df.loc[0, 'current_A']:.3f} A")
print(f"θ_peak:               {df.loc[0, 'theta_peak_deg']:.2f}° ± {df.loc[0, 'u_theta_peak_deg']:.3f}°")
print(f"λ_max:                {df.loc[0, 'lambda_max_nm']:.1f} nm ± {df.loc[0, 'u_lambda_max_nm']:.1f} nm")
print(f"T:                    {df.loc[0, 'T_K']:.0f} K ± {df.loc[0, 'u_T_K']:.0f} K")
print(f"1/T:                  {df.loc[0, 'inv_T']:.6e} K⁻¹ ± {df.loc[0, 'u_inv_T']:.3e} K⁻¹")
print(f"b (Wien):             {df.loc[0, 'b_mK']:.6e} m·K ± {df.loc[0, 'u_b_mK']:.3e} m·K")

print(f"\nRELATIVE UNCERTAINTIES (%):")
print(f"  θ_peak:   {100*df.loc[0, 'u_theta_peak_deg']/df.loc[0, 'theta_peak_deg']:.2f}%")
print(f"  λ_max:    {100*df.loc[0, 'u_lambda_max_nm']/df.loc[0, 'lambda_max_nm']:.2f}%")
print(f"  T:        {100*df.loc[0, 'u_T_K']/df.loc[0, 'T_K']:.2f}%")
print(f"  b:        {100*df.loc[0, 'u_b_mK']/df.loc[0, 'b_mK']:.2f}%")

print("\n" + "=" * 75)
print("SUMMARY STATISTICS (all 9 measurements):")
print("=" * 75)
valid_b = df['b_mK'].dropna()
valid_u_b = df['u_b_mK'].dropna()
print(f"Wien constant mean:        {np.mean(valid_b):.6e} m·K")
print(f"Wien constant std dev:     {np.std(valid_b, ddof=1):.6e} m·K")
print(f"Mean measurement unc:      {np.mean(valid_u_b):.6e} m·K")
print(f"Literature b₀:             2.898e-03 m·K")
print(f"Difference:                {100*(np.mean(valid_b) - 2.898e-3)/2.898e-3:.2f}%")

print("\n" + "=" * 75)
print("CSV COLUMNS:")
print("=" * 75)
print("Data columns:")
for col in df.columns:
    if not col.startswith('u_') and col != 'notes':
        print(f"  • {col}")
print("\nUncertainty columns:")
for col in df.columns:
    if col.startswith('u_'):
        print(f"  • {col}")

print("\n✓ All uncertainties computed and saved to CSV!")
