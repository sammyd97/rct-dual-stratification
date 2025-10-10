#!/usr/bin/env python3
"""
Randomize facilities for an RCT without stratification.
Produces reproducible allocation and balance diagnostics.
"""

import argparse
import sys
import pandas as pd
import numpy as np
from pathlib import Path

# Try importing scipy for statistical tests
try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("Warning: scipy not available. Statistical tests will be skipped.", file=sys.stderr)


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Randomize facilities for RCT without stratification"
    )
    parser.add_argument(
        "--in", dest="input_csv", required=True,
        help="Input CSV file path (UTF-8)"
    )
    parser.add_argument(
        "--out", dest="output_csv", required=True,
        help="Output CSV file path with Arm and RandID columns"
    )
    parser.add_argument(
        "--seed", type=int, default=12345,
        help="Random seed for reproducibility (default: 12345)"
    )
    parser.add_argument(
        "--arms", type=str, default="Control,Treatment",
        help="Comma-separated arm labels (default: Control,Treatment)"
    )
    parser.add_argument(
        "--ratio", type=str, default="1,1",
        help="Comma-separated allocation ratio (default: 1,1)"
    )
    parser.add_argument(
        "--report-cols", type=str, default=None,
        help="Comma-separated column names for balance diagnostics"
    )
    parser.add_argument(
        "--report-out", type=str, default=None,
        help="Optional: path to save detailed balance report CSV"
    )
    parser.add_argument(
        "--stratify-by", type=str, default=None,
        help="Optional: column name to stratify randomization by (e.g., LGA)"
    )
    return parser.parse_args()


def validate_inputs(arms, ratio, df, report_cols, stratify_by):
    """Validate input parameters."""
    # Check arms and ratio lengths match
    if len(arms) != len(ratio):
        raise ValueError(
            f"Number of arms ({len(arms)}) must match number of ratio values ({len(ratio)})"
        )
    
    # Check at least 2 arms
    if len(arms) < 2:
        raise ValueError("Must specify at least 2 arms")
    
    # Check ratio values are positive
    if any(r <= 0 for r in ratio):
        raise ValueError("All ratio values must be positive")
    
    # Check dataframe not empty
    if df.empty:
        raise ValueError("Input CSV is empty")
    
    # Check report columns exist
    if report_cols:
        missing_cols = set(report_cols) - set(df.columns)
        if missing_cols:
            raise ValueError(f"Report columns not found in data: {missing_cols}")
    
    # Check stratification column exists
    if stratify_by and stratify_by not in df.columns:
        raise ValueError(f"Stratification column '{stratify_by}' not found in data")


def perform_randomization(df, arms, ratio, seed, stratify_by=None):
    """
    Perform permuted block randomization with optional stratification.
    
    Args:
        df: Input dataframe
        arms: List of arm names
        ratio: List of allocation ratios
        seed: Random seed
        stratify_by: Optional column name to stratify by
    
    Returns:
        DataFrame with added RandID and Arm columns
    """
    rng = np.random.default_rng(seed)
    
    # Initialize output dataframe
    df_output = df.copy()
    df_output['RandID'] = 0.0
    df_output['Arm'] = ''
    
    # If no stratification, randomize the whole dataset
    if not stratify_by:
        return _randomize_group(df, df.index, arms, ratio, rng)
    
    # Stratified randomization
    print(f"\n=== Stratified Randomization by {stratify_by} ===")
    
    strata = df[stratify_by].unique()
    for stratum in sorted(strata):
        stratum_indices = df[df[stratify_by] == stratum].index
        stratum_size = len(stratum_indices)
        
        # Randomize within this stratum
        df_stratum = _randomize_group(df.loc[stratum_indices], stratum_indices, arms, ratio, rng)
        
        # Copy results to output
        df_output.loc[stratum_indices, 'RandID'] = df_stratum.loc[stratum_indices, 'RandID']
        df_output.loc[stratum_indices, 'Arm'] = df_stratum.loc[stratum_indices, 'Arm']
        
        # Print stratum allocation
        stratum_counts = df_output.loc[stratum_indices, 'Arm'].value_counts()
        print(f"{stratum} (n={stratum_size}):", end=" ")
        for arm in arms:
            count = stratum_counts.get(arm, 0)
            print(f"{arm}={count}", end=" ")
        print()
    
    print()
    
    return df_output


def _randomize_group(df, indices, arms, ratio, rng):
    """
    Randomize a single group using permuted block randomization.
    
    Args:
        df: Full dataframe
        indices: Indices of rows to randomize
        arms: List of arm names
        ratio: List of allocation ratios
        rng: Random number generator
    
    Returns:
        DataFrame with RandID and Arm assigned for the specified indices
    """
    n = len(indices)
    
    # Create RandID for each row
    rand_ids = rng.uniform(0, 1, size=n)
    
    # Create temporary dataframe with indices and RandIDs
    temp_df = pd.DataFrame({
        'orig_idx': indices,
        'RandID': rand_ids
    })
    
    # Sort by RandID for randomization
    temp_df = temp_df.sort_values('RandID').reset_index(drop=True)
    
    # Permuted block allocation
    block_size = sum(ratio)
    n_full_blocks = n // block_size
    remainder = n % block_size
    
    # Create allocation sequence
    allocation = []
    
    # Full blocks
    block_template = []
    for arm, count in zip(arms, ratio):
        block_template.extend([arm] * count)
    
    for _ in range(n_full_blocks):
        # Permute each block
        block = block_template.copy()
        rng.shuffle(block)
        allocation.extend(block)
    
    # Handle remainder approximately in proportion
    if remainder > 0:
        # Create a partial block with approximate proportions
        remainder_counts = []
        total_ratio = sum(ratio)
        allocated_so_far = 0
        
        for i, (arm, r) in enumerate(zip(arms, ratio)):
            if i < len(arms) - 1:
                # Allocate proportionally, rounding
                count = round(remainder * r / total_ratio)
                remainder_counts.append(count)
                allocated_so_far += count
            else:
                # Last arm gets whatever is left
                remainder_counts.append(remainder - allocated_so_far)
        
        remainder_block = []
        for arm, count in zip(arms, remainder_counts):
            remainder_block.extend([arm] * count)
        
        rng.shuffle(remainder_block)
        allocation.extend(remainder_block)
    
    # Assign arms to temp dataframe
    temp_df['Arm'] = allocation
    
    # Round RandID to 6 decimals
    temp_df['RandID'] = temp_df['RandID'].round(6)
    
    # Create output with original indices
    df_output = df.copy()
    df_output['RandID'] = 0.0
    df_output['Arm'] = ''
    
    for _, row in temp_df.iterrows():
        orig_idx = row['orig_idx']
        df_output.loc[orig_idx, 'RandID'] = row['RandID']
        df_output.loc[orig_idx, 'Arm'] = row['Arm']
    
    return df_output


def compute_smd(mean1, sd1, n1, mean2, sd2, n2):
    """Compute Standardized Mean Difference (Cohen's d with pooled SD)."""
    if sd1 == 0 and sd2 == 0:
        return 0.0
    pooled_sd = np.sqrt(((n1 - 1) * sd1**2 + (n2 - 1) * sd2**2) / (n1 + n2 - 2))
    if pooled_sd == 0:
        return 0.0
    return (mean1 - mean2) / pooled_sd


def generate_balance_report(df, arms, report_cols):
    """
    Generate balance diagnostics for specified columns.
    
    Returns:
        Tuple of (summary_text, detailed_dataframe)
    """
    summary_lines = []
    detail_rows = []
    
    # Overall allocation counts
    summary_lines.append("=== Allocation (overall) ===")
    summary_lines.append(f"{'Arm':<20} {'N':>8}")
    for arm in arms:
        count = (df['Arm'] == arm).sum()
        summary_lines.append(f"{arm:<20} {count:>8}")
    summary_lines.append("")
    
    if not report_cols:
        return "\n".join(summary_lines), pd.DataFrame()
    
    # Separate numeric and categorical variables
    numeric_vars = []
    categorical_vars = []
    
    for col in report_cols:
        if pd.api.types.is_numeric_dtype(df[col]):
            numeric_vars.append(col)
        else:
            categorical_vars.append(col)
    
    # Numeric balance
    if numeric_vars:
        summary_lines.append("=== Balance (numeric) ===")
        header = f"{'Variable':<30} {'Arm':<15} {'N':>6} {'Mean':>10} {'SD':>10} {'SMD':>8} {'p-value':>10}"
        summary_lines.append(header)
        
        for var in numeric_vars:
            # Collect stats for each arm
            stats_by_arm = {}
            for arm in arms:
                arm_data = df[df['Arm'] == arm][var].dropna()
                stats_by_arm[arm] = {
                    'n': len(arm_data),
                    'mean': arm_data.mean() if len(arm_data) > 0 else np.nan,
                    'sd': arm_data.std(ddof=1) if len(arm_data) > 1 else 0.0
                }
            
            # Compute SMD and p-value (for first two arms)
            if len(arms) >= 2:
                arm1, arm2 = arms[0], arms[1]
                s1, s2 = stats_by_arm[arm1], stats_by_arm[arm2]
                
                if s1['n'] > 0 and s2['n'] > 0:
                    smd = compute_smd(s1['mean'], s1['sd'], s1['n'], 
                                     s2['mean'], s2['sd'], s2['n'])
                    
                    # Welch's t-test
                    p_val = np.nan
                    if SCIPY_AVAILABLE:
                        data1 = df[df['Arm'] == arm1][var].dropna()
                        data2 = df[df['Arm'] == arm2][var].dropna()
                        if len(data1) > 0 and len(data2) > 0:
                            _, p_val = stats.ttest_ind(data1, data2, equal_var=False)
                else:
                    smd = np.nan
                    p_val = np.nan
            else:
                smd = np.nan
                p_val = np.nan
            
            # Print first arm with SMD and p-value
            arm = arms[0]
            s = stats_by_arm[arm]
            mean_str = f"{s['mean']:.1f}" if not np.isnan(s['mean']) else "NA"
            sd_str = f"{s['sd']:.1f}" if not np.isnan(s['sd']) else "NA"
            smd_str = f"{smd:.3f}" if not np.isnan(smd) else "NA"
            p_str = f"{p_val:.3f}" if not np.isnan(p_val) else "NA"
            
            summary_lines.append(
                f"{var:<30} {arm:<15} {s['n']:>6} {mean_str:>10} {sd_str:>10} {smd_str:>8} {p_str:>10}"
            )
            
            # Print other arms without SMD/p-value
            for arm in arms[1:]:
                s = stats_by_arm[arm]
                mean_str = f"{s['mean']:.1f}" if not np.isnan(s['mean']) else "NA"
                sd_str = f"{s['sd']:.1f}" if not np.isnan(s['sd']) else "NA"
                summary_lines.append(
                    f"{'':<30} {arm:<15} {s['n']:>6} {mean_str:>10} {sd_str:>10}"
                )
            
            # Add to detailed report
            for arm in arms:
                s = stats_by_arm[arm]
                detail_rows.append({
                    'variable': var,
                    'type': 'numeric',
                    'level': np.nan,
                    'arm': arm,
                    'n': s['n'],
                    'mean': s['mean'],
                    'sd': s['sd'],
                    'prop': np.nan,
                    'smd': smd if arm == arms[0] else np.nan,
                    'p_value': p_val if arm == arms[0] else np.nan
                })
        
        summary_lines.append("")
    
    # Categorical balance
    if categorical_vars:
        summary_lines.append("=== Balance (categorical) ===")
        
        for var in categorical_vars:
            # Handle missing values
            df_var = df[[var, 'Arm']].copy()
            df_var[var] = df_var[var].fillna('__MISSING__')
            
            # Get unique levels
            levels = sorted(df_var[var].unique())
            
            # Create contingency table
            ct = pd.crosstab(df_var[var], df_var['Arm'])
            
            # Chi-square test
            p_val = np.nan
            if SCIPY_AVAILABLE and len(levels) > 1 and len(arms) > 1:
                try:
                    _, p_val, _, _ = stats.chi2_contingency(ct)
                except:
                    pass
            
            # Print header for this variable
            arm_headers = []
            for arm in arms:
                arm_headers.append(f"{arm}_n")
                arm_headers.append(f"{arm}_%")
            
            header_parts = [f"{'Variable':<20}", f"{'Level':<20}"]
            for arm in arms:
                header_parts.append(f"{arm+'_n':>10}")
                header_parts.append(f"{arm+'_%':>10}")
            header_parts.append(f"{'p-value':>10}")
            header_parts.append(f"{'max|Δ%|':>10}")
            
            summary_lines.append(" ".join(header_parts))
            
            # Calculate max absolute difference in proportions
            max_diff = 0.0
            for level in levels:
                proportions = []
                for arm in arms:
                    arm_total = (df_var['Arm'] == arm).sum()
                    level_count = ((df_var[var] == level) & (df_var['Arm'] == arm)).sum()
                    prop = 100 * level_count / arm_total if arm_total > 0 else 0
                    proportions.append(prop)
                
                if len(proportions) >= 2:
                    diff = max(proportions) - min(proportions)
                    max_diff = max(max_diff, diff)
            
            # Print each level
            for i, level in enumerate(levels):
                row_parts = []
                if i == 0:
                    row_parts.append(f"{var:<20}")
                else:
                    row_parts.append(f"{'':<20}")
                
                level_display = level if level != '__MISSING__' else 'NA'
                row_parts.append(f"{level_display:<20}")
                
                for arm in arms:
                    arm_total = (df_var['Arm'] == arm).sum()
                    level_count = ((df_var[var] == level) & (df_var['Arm'] == arm)).sum()
                    prop = 100 * level_count / arm_total if arm_total > 0 else 0
                    row_parts.append(f"{level_count:>10}")
                    row_parts.append(f"{prop:>9.1f}%")
                
                # Add p-value and max diff only on first row
                if i == 0:
                    p_str = f"{p_val:.3f}" if not np.isnan(p_val) else "NA"
                    row_parts.append(f"{p_str:>10}")
                    row_parts.append(f"{max_diff:>9.1f}%")
                else:
                    row_parts.append(f"{'':>10}")
                    row_parts.append(f"{'':>10}")
                
                summary_lines.append(" ".join(row_parts))
                
                # Add to detailed report
                for arm in arms:
                    arm_total = (df_var['Arm'] == arm).sum()
                    level_count = ((df_var[var] == level) & (df_var['Arm'] == arm)).sum()
                    prop = level_count / arm_total if arm_total > 0 else 0
                    
                    detail_rows.append({
                        'variable': var,
                        'type': 'categorical',
                        'level': level_display,
                        'arm': arm,
                        'n': level_count,
                        'mean': np.nan,
                        'sd': np.nan,
                        'prop': prop,
                        'smd': np.nan,
                        'p_value': p_val if i == 0 else np.nan
                    })
            
            summary_lines.append("")
    
    detail_df = pd.DataFrame(detail_rows) if detail_rows else pd.DataFrame()
    
    return "\n".join(summary_lines), detail_df


def main():
    """Main execution function."""
    args = parse_arguments()
    
    # Parse arms and ratio
    arms = [a.strip() for a in args.arms.split(',')]
    ratio = [int(r.strip()) for r in args.ratio.split(',')]
    
    # Parse report columns
    report_cols = None
    if args.report_cols:
        report_cols = [c.strip() for c in args.report_cols.split(',')]
    
    # Read input CSV
    try:
        df = pd.read_csv(args.input_csv, encoding='utf-8')
    except Exception as e:
        print(f"Error reading input CSV: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Validate inputs
    try:
        validate_inputs(arms, ratio, df, report_cols, args.stratify_by)
    except ValueError as e:
        print(f"Validation error: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Perform randomization
    df_randomized = perform_randomization(df, arms, ratio, args.seed, args.stratify_by)
    
    # Save output
    try:
        df_randomized.to_csv(args.output_csv, index=False, encoding='utf-8')
        print(f"✓ Randomized data saved to: {args.output_csv}")
    except Exception as e:
        print(f"Error writing output CSV: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Generate and print balance report
    summary_text, detail_df = generate_balance_report(df_randomized, arms, report_cols)
    print("\n" + summary_text)
    
    # Save detailed report if requested
    if args.report_out and not detail_df.empty:
        try:
            detail_df.to_csv(args.report_out, index=False, encoding='utf-8')
            print(f"✓ Detailed balance report saved to: {args.report_out}")
        except Exception as e:
            print(f"Error writing report CSV: {e}", file=sys.stderr)
    
    print(f"\n✓ Randomization complete with seed={args.seed}")
    print(f"  Same seed will produce identical results for reproducibility.")


if __name__ == "__main__":
    main()

