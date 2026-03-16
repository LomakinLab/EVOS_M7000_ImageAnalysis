import pandas as pd
import glob
import os
import numpy as np

MIN_POSITIVE_VALUE = 1e-12

PARAMETER_MAP = {
    'Area': 'Area', 'Mean': 'MGV', 'StdDev': 'StdDev',
    'IntDen': 'ID', 'RawIntDen': 'RawIntDen', 'Perim.': 'Perimeter'
}

# Specification Lists
BASE_PARAMETERS_FOR_RATIO = ['Area', 'Mean', 'StdDev', 'IntDen', 'RawIntDen']
BASE_PARAMETERS_ALL = ['Area', 'Mean', 'StdDev', 'IntDen', 'RawIntDen', 'Perim.']
INTENSITY_PARAMETERS_RAW = ['Mean', 'IntDen', 'RawIntDen']
INTENSITY_PARAMETERS_DERIVED = [
    'MGV_nuclei', 'ID_nuclei', 'RawIntDen_nuclei', 'MGV_cytoplasm', 'ID_cytoplasm',
    'MGV_average_nuclei', 'ID_average_nuclei', 'RawIntDen_average_nuclei'
]

def safe_ratio(numerator, denominator):
    if denominator == 0 or np.isnan(numerator) or np.isnan(denominator):
        return np.nan
    return numerator / denominator

def calculate_eop(P, A):
    if P == 0 or A <= 0 or np.isnan(P) or np.isnan(A):
        return np.nan
    P_circle = 2 * np.sqrt(np.pi * A)
    return (P - P_circle) / P_circle if P_circle != 0 else np.nan

def apply_floor_correction(data_list):
    return [val if (pd.isna(val) or val >= MIN_POSITIVE_VALUE) else MIN_POSITIVE_VALUE for val in data_list]

# --- 1. GLOBAL REPLICATE CONFIGURATION ---
search_root = os.path.abspath('.')
global_cond_file = os.path.join(search_root, 'conditions.txt')
condition_to_replicates = {}

if os.path.exists(global_cond_file):
    with open(global_cond_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or ':' not in line: continue
            c_name, content = line.split(':', 1)
            clean_content = content.replace(';', '').strip()
            replicates_list = [r.strip() for r in clean_content.split(',') if r.strip()]
            condition_to_replicates[c_name.strip()] = replicates_list
    print(f"Loaded global conditions from: {global_cond_file}")
else:
    print("CRITICAL: conditions.txt not found in root.")

all_conditions = sorted(condition_to_replicates.keys())

# --- 2. EXPERIMENT DISCOVERY ---
subfolders = [f.path for f in os.scandir(search_root) if f.is_dir()]
valid_exp_folders = [f for f in subfolders if os.path.exists(os.path.join(f, 'output_ijm'))]

for exp_folder in valid_exp_folders:
    input_dir = os.path.join(exp_folder, 'output_ijm')
    output_dir = os.path.join(exp_folder, 'combined_output')
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Processing experiment: {os.path.basename(exp_folder)}")

    derived_metrics = ["#nuclei", "area_nuclei", "MGV_nuclei", "ID_nuclei", "StdDev_nuclei", "RawIntDen_nuclei", 
                      "area_cytoplasm", "MGV_cytoplasm", "ID_cytoplasm", "StdDev_cytoplasm", "Eop", "CoV_nuclei_area"]
    
    ratio_params = []
    for par in BASE_PARAMETERS_FOR_RATIO:
        m = PARAMETER_MAP.get(par, par)
        ratio_params.append(f'{m}_nuclei_cell_ratio')
        if par != 'RawIntDen': ratio_params.append(f'{m}_nuclei_cytoplasm_ratio')

    avg_nuc_params = [f'{PARAMETER_MAP.get(p, p)}_average_nuclei' for p in BASE_PARAMETERS_ALL]
    track_params = BASE_PARAMETERS_ALL + derived_metrics + ratio_params + avg_nuc_params

    all_data = {c: {p: [] for p in track_params} for c in all_conditions}
    single_nuc_data = {p: {c: [] for c in all_conditions} for p in BASE_PARAMETERS_ALL}
    single_eop_data = {c: [] for c in all_conditions}
    field_counts = {c: [] for c in all_conditions}

    all_available_folders = os.listdir(input_dir)

    # 3. PROCESS BY CONDITION
    for condition in all_conditions:
        replicate_prefixes = condition_to_replicates[condition]
        
        for prefix in replicate_prefixes:
            replicate_cell_total = 0
            target_folders = [f for f in all_available_folders if f.startswith(prefix)]
            
            # Even if target_folders is empty, we must process the logic to maintain spacer sync
            for folder in sorted(target_folders):
                folder_path = os.path.join(input_dir, folder)
                cell_files = sorted(glob.glob(os.path.join(folder_path, 'cell_*.csv')))

                for cell_file in cell_files:
                    try:
                        df_cell = pd.read_csv(cell_file)
                        if df_cell.empty: continue
                        cell_data = df_cell.iloc[0].to_dict()
                        c_idx = os.path.basename(cell_file).split('_')[1].split('.')[0]
                        
                        bg_mgv = 0.0
                        bg_path = os.path.join(folder_path, f'background_for_cell_{c_idx}.csv')
                        if os.path.exists(bg_path):
                            bg_mgv = pd.read_csv(bg_path, index_col=0).iloc[0].get('Mean', 0.0)

                        cell_data['Mean'] -= bg_mgv
                        for col in ['IntDen', 'RawIntDen']:
                            if col in cell_data: cell_data[col] -= (bg_mgv * cell_data['Area'])

                        nuc_path = os.path.join(folder_path, f'nuclei_for_cell_{c_idx}.csv')
                        df_nuc = pd.read_csv(nuc_path) if os.path.exists(nuc_path) else pd.DataFrame()
                        
                        cur_derived = {p: np.nan for p in track_params}
                        num_nucs = len(df_nuc)
                        cur_derived["#nuclei"] = num_nucs

                        if not df_nuc.empty:
                            replicate_cell_total += 1
                            df_nuc['Mean'] -= bg_mgv
                            for col in ['IntDen', 'RawIntDen']:
                                if col in df_nuc.columns: df_nuc[col] -= (bg_mgv * df_nuc['Area'])

                            if num_nucs == 1: cur_derived['CoV_nuclei_area'] = 0.0
                            elif num_nucs > 1:
                                ma, sa = df_nuc['Area'].mean(), df_nuc['Area'].std()
                                cur_derived['CoV_nuclei_area'] = (sa / ma) * 100 if ma != 0 else 0.0

                            for p in BASE_PARAMETERS_ALL:
                                single_nuc_data[p][condition].extend(df_nuc[p].tolist())
                                cur_derived[f'{PARAMETER_MAP.get(p, p)}_average_nuclei'] = df_nuc[p].mean()
                            
                            for _, r in df_nuc.iterrows():
                                single_eop_data[condition].append(calculate_eop(r.get('Perim.', 0), r.get('Area', 0)))

                            cur_derived['area_nuclei'] = df_nuc['Area'].sum()
                            cur_derived['ID_nuclei'] = df_nuc['IntDen'].sum() if 'IntDen' in df_nuc.columns else 0.0
                            cur_derived['RawIntDen_nuclei'] = df_nuc['RawIntDen'].sum() if 'RawIntDen' in df_nuc.columns else 0.0
                            cur_derived['MGV_nuclei'] = safe_ratio(cur_derived['ID_nuclei'], cur_derived['area_nuclei'])
                            cur_derived['StdDev_nuclei'] = df_nuc['StdDev'].mean() if 'StdDev' in df_nuc.columns else np.nan

                        cyto_path = os.path.join(folder_path, f'cytoplasm_for_cell_{c_idx}.csv')
                        if os.path.exists(cyto_path):
                            df_cyto = pd.read_csv(cyto_path).iloc[0]
                            cur_derived['area_cytoplasm'] = df_cyto['Area']
                            cur_derived['MGV_cytoplasm'] = df_cyto['Mean'] - bg_mgv
                            cur_derived['ID_cytoplasm'] = df_cyto['IntDen'] - (bg_mgv * df_cyto['Area'])
                            cur_derived['StdDev_cytoplasm'] = df_cyto['StdDev']

                        cur_derived['Eop'] = calculate_eop(cell_data.get('Perim.', 0), cell_data.get('Area', 0))
                        for par in BASE_PARAMETERS_FOR_RATIO:
                            m = PARAMETER_MAP.get(par, par)
                            nv = cur_derived['area_nuclei'] if par == 'Area' else cur_derived.get(f'{m}_nuclei')
                            cur_derived[f'{m}_nuclei_cell_ratio'] = safe_ratio(nv, cell_data.get(par))
                            if par != 'RawIntDen':
                                cv = cur_derived['area_cytoplasm'] if par == 'Area' else cur_derived.get(f'{m}_cytoplasm')
                                cur_derived[f'{m}_nuclei_cytoplasm_ratio'] = safe_ratio(nv, cv)

                        for p in BASE_PARAMETERS_ALL: all_data[condition][p].append(cell_data.get(p, np.nan))
                        for p in track_params:
                            if p not in BASE_PARAMETERS_ALL: all_data[condition][p].append(cur_derived[p])
                    except Exception: continue

            # --- APPLY REPLICATE SPACER (Ensures every prefix has a closing spacer) ---
            field_counts[condition].extend([replicate_cell_total, np.nan])
            for p in track_params: all_data[condition][p].append(np.nan)
            for p in BASE_PARAMETERS_ALL: single_nuc_data[p][condition].append(np.nan)
            single_eop_data[condition].append(np.nan)

    # 4. SAVE OUTPUTS
    for p in INTENSITY_PARAMETERS_RAW + INTENSITY_PARAMETERS_DERIVED:
        for c in all_conditions: all_data[c][p] = apply_floor_correction(all_data[c][p])

    for p in track_params:
        df = pd.DataFrame({c: pd.Series(all_data[c][p]) for c in all_conditions})
        fname = PARAMETER_MAP.get(p, p).replace("#", "Number_")
        df.to_csv(os.path.join(output_dir, f'combined_{fname}.csv'), index=False)

    for p in BASE_PARAMETERS_ALL:
        df = pd.DataFrame({c: pd.Series(single_nuc_data[p][c]) for c in all_conditions})
        df.to_csv(os.path.join(output_dir, f'combined_{PARAMETER_MAP.get(p, p)}_single_nuclei.csv'), index=False)

    pd.DataFrame({c: pd.Series(single_eop_data[c]) for c in all_conditions}).to_csv(os.path.join(output_dir, 'combined_Eop_single_nuclei.csv'), index=False)
    pd.DataFrame({c: pd.Series(field_counts[c]) for c in all_conditions}).to_csv(os.path.join(output_dir, 'cells_per_field.csv'), index=False)
    
    m_d = {c: pd.Series(all_data[c]['Mean']) for c in all_conditions}
    s_d = {c: pd.Series(all_data[c]['StdDev']) for c in all_conditions}
    pd.DataFrame({c: s_d[c]/m_d[c] for c in all_conditions}).to_csv(os.path.join(output_dir, 'combined_CoV_cell.csv'), index=False)

print("Process finished. Alignment gremlin fixed.")
