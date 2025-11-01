import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress, spearmanr
import numpy as np
import os
import seaborn as sns
from sklearn.linear_model import RANSACRegressor
from sklearn.preprocessing import StandardScaler

DATA_DIR = './combined_output'

FILE_MAPPING = {
    'Area nuclei': 'combined_area_nuclei.csv',
    'Area cell': 'combined_Area.csv',
    'ID nuclei': 'combined_ID_nuclei.csv',
    'ID cell': 'combined_ID.csv',
    'Number nuclei': 'combined_Number_nuclei.csv',
    'Area average nuclei': 'combined_Area_average_nuclei.csv',
}

def load_and_combine_data(file_mapping, data_dir):
    all_data = {}
    all_features = file_mapping.keys()
    
    print(f"Attempting to load data from directory: {data_dir}")

    for feature_name, filename in file_mapping.items():
        filepath = os.path.join(data_dir, filename)
        
        try:
            df_wide = pd.read_csv(filepath)
            
            df_wide['Measurement_ID'] = df_wide.index
            
            df_long = pd.melt(
                df_wide,
                id_vars=['Measurement_ID'],
                var_name='Condition',
                value_name=feature_name
            )
            
            df_long = df_long.set_index(['Measurement_ID', 'Condition'])
            
            all_data[feature_name] = df_long
            
        except FileNotFoundError:
            print(f"FATAL ERROR: Required file not found: {filepath}.")
            return None
        except Exception as e:
            print(f"FATAL ERROR while processing {filepath}: {e}")
            return None

    if not all_data or len(all_data) != len(all_features):
        print("Error: Did not successfully load all required feature files.")
        return None
        
    master_df = pd.concat(all_data.values(), axis=1).reset_index()

    master_df.dropna(subset=['Condition'], inplace=True)
    
    return master_df


def calculate_bivariate_stats(data_df, x_col, y_col):
    
    sub_df = data_df[[x_col, y_col]].dropna().copy()
    
    if len(sub_df) < 3:
        return None
        
    X_full = sub_df[x_col].values.reshape(-1, 1)
    Y_full = sub_df[y_col].values
    
    original_n = len(X_full)

    try:
        scaler_x = StandardScaler()
        scaler_y = StandardScaler()
        X_scaled = scaler_x.fit_transform(X_full)
        Y_scaled = scaler_y.fit_transform(Y_full.reshape(-1, 1)).flatten()
        
        ransac = RANSACRegressor(
            min_samples=2,
            max_trials=1000,
            random_state=42,
            residual_threshold=2.0
        )
        ransac.fit(X_scaled, Y_scaled)
        
        inlier_mask = ransac.inlier_mask_
        
        X_inliers_scaled = X_scaled[inlier_mask]
        Y_inliers_scaled = Y_scaled[inlier_mask]
        
        slope_scaled = ransac.estimator_.coef_[0]
        intercept_scaled = ransac.estimator_.intercept_
        
        slope = slope_scaled * (scaler_y.scale_[0] / scaler_x.scale_[0])
        intercept = scaler_y.mean_[0] + intercept_scaled * scaler_y.scale_[0] - slope * scaler_x.mean_[0]
        
        X = X_full[inlier_mask].flatten()
        Y = Y_full[inlier_mask]
        
        line = slope * X + intercept
        
        n_points_filtered = len(X)
        n_removed = original_n - n_points_filtered

    except Exception as e:
        print(f"RANSAC failed ({e}). Falling back to standard linregress on all data.")
        slope, intercept, _, _, _ = linregress(X_full.flatten(), Y_full)
        X, Y = X_full.flatten(), Y_full
        line = slope * X + intercept
        n_points_filtered = original_n
        n_removed = 0

    try:
        rho, p_value_spearman = spearmanr(X, Y)
    except ValueError:
        rho, p_value_spearman = np.nan, np.nan

    return {
        'n_points': n_points_filtered,
        'original_n': original_n,
        'X': X,
        'Y': Y,
        'slope': slope,
        'intercept': intercept,
        'rho': rho,
        'p_value_raw': p_value_spearman,
        'line': line,
        'n_removed': n_removed
    }

def plot_correlation_per_condition(df_full, x_col, y_col, condition, global_limits):
    
    output_dir = './correlation_plots/'
    os.makedirs(output_dir, exist_ok=True)
    
    df_condition = df_full[df_full['Condition'] == condition]
    
    stats = calculate_bivariate_stats(df_condition, x_col, y_col)
    
    if stats is None or np.isnan(stats.get('rho', np.nan)):
        print(f"Skipping {y_col} vs {x_col} for {condition}: Insufficient data (n < 3) or correlation failed.")
        return

    x_lim = global_limits['x_lim']
    y_lim = global_limits['y_lim']

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    ax.scatter(stats['X'], stats['Y'], s=50, alpha=0.7, edgecolors='k', color=sns.color_palette("Set1")[0],
               label=f"Measurements (n={stats['n_points']})")
    
    ax.plot(stats['X'], stats['line'], color='red', linestyle='-', linewidth=2)
    
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)

    sig_color = 'green' if stats['p_value_raw'] < 0.05 else 'black'
    
    n_removed = stats['n_removed']

    slope_value = stats["slope"]
    if abs(slope_value) > 0.01:
        slope_formatted = f'{slope_value:.3f}'
    else:
        slope_formatted = f'{slope_value:.3e}'

    stat_text = (
        f'n (used): {stats["n_points"]} (removed: {n_removed})\n'
        f'Slope (RANSAC): {slope_formatted}\n'
        r'Spearman $\rho$: ' + f'{stats["rho"]:.3f}\n'
        r'P-value: ' + f'{stats["p_value_raw"]:.3e}'
    )
    
    ax.text(0.95, 0.05, stat_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment='bottom',
            horizontalalignment='right',
            color=sig_color,
            bbox=dict(boxstyle="round,pad=0.5", fc="white", alpha=0.9,
                      edgecolor=sig_color if stats['p_value_raw'] < 0.05 else 'gray'))

    ax.set_title(f'{y_col} vs. {x_col} ({condition}) - RANSAC Filtered', fontsize=14, fontweight='bold')
    ax.set_xlabel(x_col, fontsize=12)
    ax.set_ylabel(y_col, fontsize=12)
    
    ax.grid(True, linestyle='--', alpha=0.5)

    safe_x_col = x_col.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '_')
    safe_y_col = y_col.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '_')
    safe_condition = condition.replace(' ', '_').replace('(', '').replace(')', '').replace('/', '_')
    
    filename = f'{safe_y_col}_vs_{safe_x_col}_{safe_condition}_RANSAC_Filtered.png'
    plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved plot: {filename}")

if __name__ == '__main__':
    
    CORRELATION_PAIRS = [
        ('Area nuclei', 'Area cell'),
        ('ID nuclei', 'Area nuclei'),
        ('ID nuclei', 'Area cell'),
        ('ID cell', 'Area cell'),
        ('Number nuclei', 'Area cell'),
        ('Number nuclei', 'Area nuclei'),
        ('Number nuclei', 'Area average nuclei'),
    ]
    
    data_df = load_and_combine_data(FILE_MAPPING, DATA_DIR)

    if data_df is None:
        print("\nAnalysis terminated due to critical data loading errors. Please ensure all files exist in the correct directory.")
        exit()
        
    print(f"\nSuccessfully loaded and combined data for {len(data_df)} individual measurements.")

    required_cols = list(FILE_MAPPING.keys())
    required_cols.append('Condition')
    
    missing_cols = [col for col in required_cols if col not in data_df.columns]
    
    if missing_cols:
        print(f"\nFATAL ERROR: The following required columns are missing after merging: {', '.join(missing_cols)}")
        print("This suggests an issue with the file mapping or the data structure.")
        exit()

    unique_conditions = data_df['Condition'].unique()
    print(f"Found {len(unique_conditions)} unique conditions: {', '.join(unique_conditions)}")
    
    unique_conditions = [c for c in unique_conditions if c is not None and c == c]

    print("\nStarting Two-Pass Analysis: Pre-calculating global axis limits based on RANSAC inliers...")
    GLOBAL_INLIER_LIMITS = {}
    
    for y_col, x_col in CORRELATION_PAIRS:
        all_x_inliers = []
        all_y_inliers = []
        pair_key = (x_col, y_col)

        for condition in unique_conditions:
            df_condition = data_df[data_df['Condition'] == condition]
            stats = calculate_bivariate_stats(df_condition, x_col, y_col)
            
            if stats is not None:
                all_x_inliers.extend(stats['X'])
                all_y_inliers.extend(stats['Y'])

        if all_x_inliers and all_y_inliers:
            x_min_global = np.min(all_x_inliers)
            x_max_global = np.max(all_x_inliers)
            y_min_global = np.min(all_y_inliers)
            y_max_global = np.max(all_y_inliers)

            x_range = x_max_global - x_min_global
            y_range = y_max_global - y_min_global

            x_padding = x_range * 0.05
            y_padding = y_range * 0.05
            
            x_lim = (x_min_global - x_padding, x_max_global + x_padding)
            y_lim = (y_min_global - y_padding, y_max_global + y_padding)

            GLOBAL_INLIER_LIMITS[pair_key] = {
                'x_lim': x_lim,
                'y_lim': y_lim
            }
        else:
             print(f"Warning: No valid RANSAC inliers found for pair: {y_col} vs {x_col}. Falling back to full data range.")
             x_min_global = data_df[x_col].min()
             x_max_global = data_df[x_col].max()
             x_range = x_max_global - x_min_global
             x_padding = x_range * 0.05
             
             y_min_global = data_df[y_col].min()
             y_max_global = data_df[y_col].max()
             y_range = y_max_global - y_min_global
             y_padding = y_range * 0.05
             
             GLOBAL_INLIER_LIMITS[pair_key] = {
                'x_lim': (x_min_global - x_padding, x_max_global + x_padding),
                'y_lim': (y_min_global - y_padding, y_max_global + y_padding)
               }


    print("\n--- Starting Plotting Phase ---")
    
    for y_col, x_col in CORRELATION_PAIRS:
        pair_key = (x_col, y_col)
        limits = GLOBAL_INLIER_LIMITS.get(pair_key)
        
        if limits is None:
             print(f"Skipping plotting for {y_col} vs {x_col}: Limits could not be determined.")
             continue
            
        for condition in unique_conditions:
            plot_correlation_per_condition(data_df, x_col, y_col, condition, limits)

    print("\n--- All correlation analysis and plotting complete ---")
    print("Plots saved to the 'correlation_plots' folder.")