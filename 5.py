import pandas as pd
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests
from itertools import combinations
import numpy as np

plot_significance = True
AGGRESSIVE_CAP_UPPER_Q = 0.95
AGGRESSIVE_CAP_LOWER_Q = 0.05
NORMAL_CAP_UPPER_Q = 0.99
NORMAL_CAP_LOWER_Q = 0.01

def add_significance_stars(ax, pairs, p_values, conditions_list, y_start, h, col):
    """
    Adds significance lines and stars for multiple pairwise comparisons (e.g., Tukey, Bonferroni-corrected M-W).
    """
    y_offset = 0
    y_increment = 6.0 * h
    
    x_coords = {cond: i for i, cond in enumerate(conditions_list)}
    
    significant_pairs = sorted([(pair, p) for pair, p in zip(pairs, p_values) if p < 0.05],
                               key=lambda x: (abs(x_coords[x[0][1]] - x_coords[x[0][0]]),
                                              min(x_coords[x[0][0]], x_coords[x[0][1]])))

    sig_count = 0
    
    for pair, p_val in significant_pairs:
        cond1, cond2 = pair
        
        if cond1 in x_coords and cond2 in x_coords:
            x1 = x_coords[cond1]
            x2 = x_coords[cond2]
            
            if x1 > x2:
                x1, x2 = x2, x1
            
            if p_val < 0.001:
                label = '***'
            elif p_val < 0.01:
                label = '**'
            else:
                label = '*'
            
            y = y_start + y_offset
            ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c='k')
            
            ax.text((x1 + x2) * .5, y + h, label, ha='center', va='bottom', color='k', fontsize=12)
            
            y_offset += y_increment
            sig_count += 1
            
    return sig_count, y_offset


def analyze_and_plot(df_param, param_name, output_dir):
    
    raw_long_df = df_param.melt(var_name='Condition', value_name='Value').dropna(subset=['Value'])
    
    is_ratio = "ratio" in param_name.lower()
    is_cov = "cov" in param_name.lower()
    
    analysis_df = raw_long_df.copy()
    plot_df = raw_long_df.copy()
    
    long_df = plot_df
    stats_df = analysis_df
    
    conditions_list = df_param.columns.tolist()
    n_conditions = len(conditions_list)
    
    if long_df.empty or n_conditions < 2:
        print(f"Skipping plot: Insufficient valid data or conditions ({n_conditions}).")
        return
    
    colors = sns.color_palette("viridis", n_colors=n_conditions)

    is_normal = True
    for cond in conditions_list:
        data = stats_df[stats_df['Condition'] == cond]['Value']
        if len(data) >= 3 and stats.shapiro(data)[1] < 0.05:
            is_normal = False
            break
        
    print(f"Normality assumption check (on analysis data): {'Passed' if is_normal else 'Failed'} (for groups with N>=3)")

    
    if is_ratio or is_cov:
        cap_upper_q = AGGRESSIVE_CAP_UPPER_Q
        cap_lower_q = AGGRESSIVE_CAP_LOWER_Q
    elif not is_normal:
        cap_upper_q = NORMAL_CAP_UPPER_Q
        cap_lower_q = NORMAL_CAP_LOWER_Q
    else:
        cap_upper_q = 1.0
        cap_lower_q = 0.0
        
    q_upper = long_df['Value'].quantile(cap_upper_q)
    q_lower = long_df['Value'].quantile(cap_lower_q)

    y_true_min = long_df['Value'].min()
    
    y_max_capped = q_upper
    y_min_fit = q_lower
    
    y_min_fit = min(y_min_fit, y_true_min)

    cap_msg = f"Using {int(cap_lower_q*100)}th to {int(cap_upper_q*100)}th percentile range for CAPPING Y-AXIS."
    print(f"Visualization Capping: {cap_msg}")
    
    y_range_capped = y_max_capped - y_min_fit
    h = y_range_capped * 0.01
    
    y_start = y_max_capped + (h * 2.5)
    
    p_val = np.nan
    ax = None
    sig_count = 0
    final_y_offset = 0

    
    try:
        plt.figure(figsize=(n_conditions * 3, 6))
        
        if is_normal:
            plot_type = 'Barplot'
            ax = sns.barplot(x='Condition', y='Value', data=long_df, palette=colors, capsize=0.1, errorbar='se', order=conditions_list)
            print("Plotting: Barplot (Normal Data)")
        else:
            plot_type = 'Boxenplot'
            ax = sns.boxenplot(x='Condition', y='Value', data=long_df, palette=colors, order=conditions_list, showfliers=False)
            print("Plotting: Boxenplot (Non-Normal Data)")
        
        if n_conditions == 2:
            data1 = stats_df[stats_df['Condition'] == conditions_list[0]]['Value']
            data2 = stats_df[stats_df['Condition'] == conditions_list[1]]['Value']
            
            if is_normal:
                t_stat, p_val = stats.ttest_ind(data1, data2, equal_var=False, nan_policy='omit')
            else:
                u_stat, p_val = stats.mannwhitneyu(data1, data2, alternative='two-sided', nan_policy='omit')
                
            if p_val < 0.05 and plot_significance:
                y_max_needed = y_start + h + (h * 0.1)
                
                if p_val < 0.001: label = '***'
                elif p_val < 0.01: label = '**'
                else: label = '*'
                
                x1, x2 = 0, 1
                ax.plot([x1, x1, x2, x2], [y_start, y_start + h, y_start + h, y_start], lw=1.5, c='k')
                ax.text((x1 + x2) * .5, y_start + h, label, ha='center', va='bottom', color='k', fontsize=12)
                sig_count = 1
            else:
                y_max_needed = y_max_capped

        elif n_conditions > 2:
            if is_normal:
                model = ols('Value ~ C(Condition)', data=stats_df).fit()
                anova_results = sm.stats.anova_lm(model, typ=2)
                
                if plot_significance and anova_results.iloc[0]['PR(>F)'] < 0.05:
                    tukey_results = pairwise_tukeyhsd(endog=stats_df['Value'], groups=stats_df['Condition'], alpha=0.05)
                    tukey_df = pd.DataFrame(data=tukey_results._results_table.data[1:], columns=tukey_results._results_table.data[0])
                    
                    significant_tukey = tukey_df[tukey_df['reject'] == True]
                    
                    if not significant_tukey.empty:
                        p_values = significant_tukey['p-adj'].values
                        pairs = list(zip(significant_tukey['group1'], significant_tukey['group2']))
                        sig_count, final_y_offset = add_significance_stars(ax, pairs, p_values, conditions_list, y_start, h, 'k')
                
            else:
                kruskal_result = stats.kruskal(*[stats_df[stats_df['Condition'] == cond]['Value'] for cond in conditions_list], nan_policy='omit')
                
                if kruskal_result.pvalue < 0.05 and plot_significance:
                    print("Kruskal-Wallis significant. Performing Pairwise Mann-Whitney U tests with Bonferroni correction.")
                    
                    p_values = []
                    pairs = list(combinations(conditions_list, 2))
                    
                    for cond1, cond2 in pairs:
                        data1 = stats_df[stats_df['Condition'] == cond1]['Value']
                        data2 = stats_df[stats_df['Condition'] == cond2]['Value']
                        u_stat, p_val = stats.mannwhitneyu(data1, data2, alternative='two-sided', nan_policy='omit')
                        p_values.append(p_val)
                        
                    reject, p_adjusted, _, _ = multipletests(p_values, alpha=0.05, method='bonferroni')
                    
                    significant_pairs_indices = np.where(reject)[0]
                    significant_pairs = [pairs[i] for i in significant_pairs_indices]
                    significant_p_values = [p_adjusted[i] for i in significant_pairs_indices]
                    
                    if significant_pairs:
                        sig_count, final_y_offset = add_significance_stars(ax, significant_pairs, significant_p_values, conditions_list, y_start, h, 'k')
                    else:
                        print("No pairwise comparisons were significant after Bonferroni correction.")
                
            y_max_needed = y_start + final_y_offset + h
        
        else:
            return
    
    except Exception as e:
        print(f"\nERROR during analysis/plotting of {param_name}: {e}")
        if ax:
            plt.close(ax.figure)
        return


    y_visual_max = y_max_capped
    y_visual_min = y_min_fit
    
    visible_range = y_visual_max - y_visual_min
    proportional_padding = visible_range * 0.05
    
    if plot_significance and sig_count > 0:
        final_top_limit = max(y_visual_max, y_max_needed) + proportional_padding
        
        ax.set_ylim(bottom=y_visual_min - proportional_padding, top=final_top_limit)

    else:
        ax.set_ylim(bottom=y_visual_min - proportional_padding, top=y_visual_max + proportional_padding)
        
    ax.set_title(f'{param_name} ({plot_type})', fontsize=14)
    ax.set_ylabel(f'{param_name}')
    ax.set_xlabel("Condition")
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{param_name}_{plot_type.lower()}.png'))
    plt.close()
    
csv_folder = './combined_output/'
plot_output_dir = './analysis_plots/'
os.makedirs(plot_output_dir, exist_ok=True)

csv_files = glob.glob(csv_folder + '*.csv')

if not csv_files:
    print(f"No CSV files found in the folder: {csv_folder}")
else:
    for file in csv_files:
        try:
            if 'combined_' in os.path.basename(file):
                    param_name = os.path.basename(file).split('combined_')[1].replace('.csv', '').replace('_', ' ')
            else:
                    param_name = os.path.basename(file).replace('.csv', '').replace('_', ' ')
                    print(f"Warning: File '{os.path.basename(file)}' does not match the 'combined_' naming convention. Using '{param_name}' as parameter name.")

        except IndexError:
            param_name = os.path.basename(file).replace('.csv', '').replace('_', ' ')
            print(f"Warning: Could not extract name properly. Using '{param_name}' as parameter name.")

        df_param = pd.read_csv(file)
        df_param.columns = df_param.columns.str.strip()
            
        print(f"\n--- Analyzing {param_name} ---")
        analyze_and_plot(df_param, param_name, plot_output_dir)

    print("\nAnalysis and plotting complete. Plots saved to the 'analysis_plots' folder.")