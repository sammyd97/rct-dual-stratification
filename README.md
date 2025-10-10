# RCT Facility Randomization Tool

This Python script performs reproducible randomization of facilities for Randomized Controlled Trials (RCTs) without stratification, and generates balance diagnostics.

## Quick Start

### Installation

You'll need Python 3 with these packages:
```bash
pip install pandas numpy scipy
```

### Basic Usage

```bash
python randomize_facilities.py \
  --in "Benue facility selection August 2025 - Copy of Sheet2.csv" \
  --out "Benue_randomised.csv" \
  --seed 20251010 \
  --arms "Control,Treatment" \
  --ratio "1,1" \
  --report-cols "LGA,2024 BCG volume" \
  --report-out "Benue_balance.csv"
```

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--in` | Yes | - | Input CSV file path |
| `--out` | Yes | - | Output CSV file path |
| `--seed` | No | 12345 | Random seed for reproducibility |
| `--arms` | No | "Control,Treatment" | Comma-separated arm labels |
| `--ratio` | No | "1,1" | Allocation ratio (must match number of arms) |
| `--report-cols` | No | None | Columns to include in balance report |
| `--report-out` | No | None | Path to save detailed balance CSV |

## Output

### 1. Randomized CSV (`--out`)
- Contains all original columns
- Plus `RandID`: random number (0-1) used for sorting
- Plus `Arm`: assigned treatment arm
- Maintains original row order

### 2. Console Output
Shows:
- Overall allocation counts by arm
- Balance for numeric variables (mean, SD, SMD, p-value)
- Balance for categorical variables (counts, proportions, chi-square test)

### 3. Detailed Balance Report (`--report-out`)
Long-format CSV with columns:
- `variable`: variable name
- `type`: "numeric" or "categorical"
- `level`: for categorical variables
- `arm`: treatment arm
- `n`: count
- `mean`, `sd`: for numeric variables
- `prop`: proportion for categorical levels
- `smd`: Standardized Mean Difference
- `p_value`: statistical test p-value

## How It Works

1. **Randomization Method**: Permuted block randomization
   - Block size = sum of allocation ratios
   - Each block is randomly permuted
   - Ensures balanced allocation

2. **Reproducibility**: Same seed → identical results
   - Change the seed to get a different (but equally valid) randomization
   - Document your seed for transparency

3. **Balance Assessment**:
   - **Numeric variables**: Welch's t-test, Standardized Mean Difference
   - **Categorical variables**: Chi-square test, proportion differences

## Examples

### 1:1 Randomization (50 facilities to Control vs Treatment)
```bash
python randomize_facilities.py \
  --in facilities.csv \
  --out randomized.csv \
  --seed 20251010 \
  --arms "Control,Treatment" \
  --ratio "1,1"
```

### 2:1 Randomization (2 Treatment for every 1 Control)
```bash
python randomize_facilities.py \
  --in facilities.csv \
  --out randomized.csv \
  --seed 20251010 \
  --arms "Control,Treatment" \
  --ratio "1,2"
```

### Three-arm Trial (1:1:1)
```bash
python randomize_facilities.py \
  --in facilities.csv \
  --out randomized.csv \
  --seed 20251010 \
  --arms "Control,Treatment A,Treatment B" \
  --ratio "1,1,1"
```

## Understanding the Seed

A **seed** is a starting number for the random number generator. Think of it as:
- **Same seed** = Same randomization (reproducible)
- **Different seed** = Different randomization (but equally valid)
- Both are unbiased and scientifically sound

**Best Practice**: Document your seed in your protocol/analysis plan.

## Tips

1. **Check your balance**: Review the balance report to ensure no major imbalances
2. **Test with different seeds**: If you see large imbalances, try a few different seeds
3. **Document everything**: Save the seed, command used, and date of randomization
4. **Run once**: Don't re-randomize unless there's a good reason (document why)

## Requirements

- Python 3.6+
- pandas
- numpy
- scipy (optional but recommended for statistical tests)

If scipy is not available, the script will still work but skip p-values.

## Support

For issues or questions, check that:
- Your input CSV has headers
- Column names in `--report-cols` match exactly (case-sensitive)
- Number of arms matches number of ratio values

---

**Example Output from Benue Data (60 facilities, seed=20251010):**

```
=== Allocation (overall) ===
Arm                         N
Control                    30
Treatment                  30

=== Balance (numeric) ===
Variable                       Arm                  N       Mean         SD      SMD    p-value
2024 BCG volume                Control             30      471.7      497.6   -0.023      0.931
                               Treatment           30      483.6      557.8

=== Balance (categorical) ===
Variable             Level                 Control_n  Control_% Treatment_n Treatment_%    p-value    max|Δ%|
LGA                  Gboko                        14      46.7%         19      63.3%      0.299      16.7%
                     Makurdi                      16      53.3%         11      36.7%
```

✅ **Good balance**: p-values > 0.05, SMD close to 0, similar proportions between arms

