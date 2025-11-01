import pandas as pd
import glob
import os
import numpy as np

PARAMETER_MAP = {
    'Area': 'Area',
    'Mean': 'MGV',
    'StdDev': 'StdDev',
    'IntDen': 'ID',
    'RawIntDen': 'RawIntDen',
    'Perim.': 'Perimeter'
}

BASE_PARAMETERS_FOR_RATIO = ['Area', 'Mean', 'StdDev', 'IntDen', 'RawIntDen']
BASE_PARAMETERS_ALL = ['Area', 'Mean', 'StdDev', 'IntDen', 'RawIntDen', 'Perim.']

input_dir = './output_ijm/'
output_dir = './combined_output/'
os.makedirs(output_dir, exist_ok=True)

folder_map = {}
conditions_file = 'conditions.txt'
if os.path.exists(conditions_file):
    with open(conditions_file, 'r') as f:
        lines = f.readlines()
        
    for line in lines:
        if ':' in line:
            condition_name, folder_prefixes_str = line.strip().split(':', 1)
            folder_prefixes = [prefix.strip() for prefix in folder_prefixes_str.replace(';', '').split(',')]
            
            for folder in os.listdir(input_dir):
                if os.path.isdir(os.path.join(input_dir, folder)):
                    for prefix in folder_prefixes:
                        if folder.upper().startswith(prefix.upper()):
                            folder_map[folder] = condition_name
                            break
else:
    print(f"Error: '{conditions_file}' not found. Please create it with the specified format.")
    folders = [f for f in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, f))]
    for folder in folders:
        if folder.upper().startswith('C'):
            folder_map[folder] = 'Ctrl'
        elif folder.upper().startswith('D'):
            folder_map[folder] = 'DTX'

all_conditions = sorted(list(set(folder_map.values())))
print(f"Detected conditions: {all_conditions}")

parameters = []

derived_parameters = [
    "#nuclei", "area_nuclei", "MGV_nuclei", "ID_nuclei", "StdDev_nuclei", "RawIntDen_nuclei",
    "area_cytoplasm", "MGV_cytoplasm", "ID_cytoplasm", "StdDev_cytoplasm",
    "Eop"
]

ratio_parameters = []
for par in BASE_PARAMETERS_FOR_RATIO:
    mapped_name = PARAMETER_MAP.get(par, par)
    ratio_parameters.append(f'{mapped_name}_nuclei_cell_ratio')
    if par != 'RawIntDen':
        ratio_parameters.append(f'{mapped_name}_nuclei_cytoplasm_ratio')

avg_nuclei_parameters = []
for par in BASE_PARAMETERS_ALL:
    mapped_name = PARAMETER_MAP.get(par, par)
    avg_nuclei_parameters.append(f'{mapped_name}_average_nuclei')

all_derived_metrics = derived_parameters + ratio_parameters + avg_nuclei_parameters

if os.path.exists(input_dir) and os.listdir(input_dir):
    first_folder = next((f for f in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, f))), None)
    if first_folder:
        first_folder_path = os.path.join(input_dir, first_folder)
        first_cell_file = glob.glob(os.path.join(first_folder_path, 'cell_*.csv'))
        
        if first_cell_file:
            df_temp = pd.read_csv(first_cell_file[0])
            if 'Unnamed: 0' in df_temp.columns:
                df_temp.drop(columns=['Unnamed: 0'], inplace=True)
            parameters = [col for col in df_temp.columns if not col.startswith('Unnamed') and col.strip()]

all_data_by_condition = {
    cond: {par: [] for par in parameters + all_derived_metrics} for cond in all_conditions
}

single_nuclei_data_by_parameter = {
    'Area': {cond: [] for cond in all_conditions},
    'Mean': {cond: [] for cond in all_conditions},
    'StdDev': {cond: [] for cond in all_conditions},
    'IntDen': {cond: [] for cond in all_conditions},
    'RawIntDen': {cond: [] for cond in all_conditions},
    'Perim.': {cond: [] for cond in all_conditions}
}

single_Eop_data_by_condition = {cond: [] for cond in all_conditions}

field_cell_counts = {}

def safe_ratio(numerator, denominator):
    if denominator == 0 or np.isnan(numerator) or np.isnan(denominator):
        return np.nan
    return numerator / denominator

def calculate_eop(P, A):
    
    if P == 0 or A <= 0 or np.isnan(P) or np.isnan(A):
        return np.nan
        
    P_circle = 2 * np.sqrt(np.pi * A)
    
    if P_circle == 0:
        return np.nan
        
    eop = (P - P_circle) / P_circle
    return eop


for folder, condition in folder_map.items():
    print(f"Processing folder: {folder} ({condition})")
    
    folder_path = os.path.join(input_dir, folder)
    cell_files = sorted(glob.glob(os.path.join(folder_path, 'cell_*.csv')))
    
    processed_cells_in_field = 0
    
    for cell_file in cell_files:
        
        derived_values = {par: np.nan if par in ratio_parameters or par == 'Eop' or par in avg_nuclei_parameters else 0.0 for par in all_derived_metrics}
        
        try:
            df_cell = pd.read_csv(cell_file)
            if 'Unnamed: 0' in df_cell.columns:
                df_cell.drop(columns=['Unnamed: 0'], inplace=True)
            cell_data = df_cell.iloc[0].to_dict()
        except pd.errors.EmptyDataError:
            print(f"Warning: Skipping empty cell file {cell_file}")
            continue

        cell_index = os.path.basename(cell_file).split('_')[1].split('.')[0]

        background_file_name = f'background_for_cell_{cell_index}.csv'
        background_file_path = os.path.join(folder_path, background_file_name)

        background_mgv = np.nan
        if os.path.exists(background_file_path):
            try:
                df_bg = pd.read_csv(background_file_path)
                if not df_bg.empty and 'Mean' in df_bg.columns:
                    background_mgv = df_bg.iloc[0].get('Mean', np.nan)
            except pd.errors.EmptyDataError:
                background_mgv = np.nan
        
        if np.isnan(background_mgv):
            background_mgv = 0.0
            print(f"Warning: Background MGV not found or empty for cell {cell_index} in {folder}. Skipping background subtraction for this cell.")


        if 'Mean' in cell_data and background_mgv != 0.0:
            cell_data['Mean'] -= background_mgv
            
            cell_area = cell_data.get('Area', 0.0)
            background_id_subtraction = background_mgv * cell_area
            
            if 'IntDen' in cell_data:
                cell_data['IntDen'] -= background_id_subtraction
            if 'RawIntDen' in cell_data:
                cell_data['RawIntDen'] -= background_id_subtraction


        if 'Perim.' in cell_data and 'Area' in cell_data:
            P_cell = cell_data['Perim.']
            A_cell = cell_data['Area']
            derived_values['Eop'] = calculate_eop(P_cell, A_cell)
        else:
            derived_values['Eop'] = np.nan
            
        nuclei_file_name = f'nuclei_for_cell_{cell_index}.csv'
        nuclei_file_path = os.path.join(folder_path, nuclei_file_name)
        
        df_nuclei = pd.DataFrame()
        if os.path.exists(nuclei_file_path):
            try:
                df_nuclei = pd.read_csv(nuclei_file_path)
                if 'Unnamed: 0' in df_nuclei.columns:
                    df_nuclei.drop(columns=['Unnamed: 0'], inplace=True)
            except pd.errors.EmptyDataError:
                df_nuclei = pd.DataFrame()


        if not df_nuclei.empty and background_mgv != 0.0:
            df_nuclei['Mean'] -= background_mgv
            
            background_id_subtraction_series = background_mgv * df_nuclei['Area']
            
            if 'IntDen' in df_nuclei.columns:
                df_nuclei['IntDen'] -= background_id_subtraction_series
            if 'RawIntDen' in df_nuclei.columns:
                df_nuclei['RawIntDen'] -= background_id_subtraction_series


        derived_values["#nuclei"] = len(df_nuclei)
        
        if derived_values["#nuclei"] > 0:
            processed_cells_in_field += 1

            for _, row in df_nuclei.iterrows():
                for par_raw in BASE_PARAMETERS_ALL:
                    if par_raw in row:
                        single_nuclei_data_by_parameter[par_raw][condition].append(row[par_raw])
                
                P_nuc = row.get('Perim.', np.nan)
                A_nuc = row.get('Area', 0.0)
                eop_nuc = calculate_eop(P_nuc, A_nuc)
                single_Eop_data_by_condition[condition].append(eop_nuc)
            
            derived_values['area_nuclei'] = df_nuclei['Area'].sum()
            derived_values['ID_nuclei'] = df_nuclei['IntDen'].sum()
            derived_values['RawIntDen_nuclei'] = df_nuclei['RawIntDen'].sum() if 'RawIntDen' in df_nuclei.columns else 0.0
            
            derived_values['MGV_nuclei'] = safe_ratio(derived_values['ID_nuclei'], derived_values['area_nuclei'])
            derived_values['StdDev_nuclei'] = df_nuclei['StdDev'].mean() if 'StdDev' in df_nuclei.columns else np.nan
            
            for par_raw in BASE_PARAMETERS_ALL:
                if par_raw in df_nuclei.columns:
                    avg_val = df_nuclei[par_raw].mean()
                    mapped_name = PARAMETER_MAP.get(par_raw, par_raw)
                    avg_param_name = f'{mapped_name}_average_nuclei'
                    derived_values[avg_param_name] = avg_val

        cytoplasm_file_name = f'cytoplasm_for_cell_{cell_index}.csv'
        cytoplasm_file_path = os.path.join(folder_path, cytoplasm_file_name)
        
        cytoplasm_data = {}
        if os.path.exists(cytoplasm_file_path):
            try:
                df_cytoplasm = pd.read_csv(cytoplasm_file_path)
                if 'Unnamed: 0' in df_cytoplasm.columns:
                    df_cytoplasm.drop(columns=['Unnamed: 0'], inplace=True)
                cytoplasm_data = df_cytoplasm.iloc[0].to_dict()
            except pd.errors.EmptyDataError:
                cytoplasm_data = {}
            
            if cytoplasm_data and background_mgv != 0.0:
                cytoplasm_data['Mean'] -= background_mgv
                
                cyto_area = cytoplasm_data.get('Area', 0.0)
                background_id_subtraction = background_mgv * cyto_area
                
                if 'IntDen' in cytoplasm_data:
                    cytoplasm_data['IntDen'] -= background_id_subtraction
                if 'RawIntDen' in cytoplasm_data:
                    cytoplasm_data['RawIntDen'] -= background_id_subtraction


            if cytoplasm_data:
                derived_values['area_cytoplasm'] = cytoplasm_data.get('Area', np.nan)
                derived_values['MGV_cytoplasm'] = cytoplasm_data.get('Mean', np.nan)
                derived_values['ID_cytoplasm'] = cytoplasm_data.get('IntDen', np.nan)
                derived_values['StdDev_cytoplasm'] = cytoplasm_data.get('StdDev', np.nan)


        for par in BASE_PARAMETERS_FOR_RATIO:
            mapped_name = PARAMETER_MAP.get(par, par)
            
            cell_component = cell_data.get(par, np.nan)
            
            if par == 'Area':
                nucleus_component = derived_values['area_nuclei']
            elif par == 'Mean':
                nucleus_component = derived_values['MGV_nuclei']
            elif par == 'StdDev':
                nucleus_component = derived_values['StdDev_nuclei']
            elif par == 'IntDen':
                nucleus_component = derived_values['ID_nuclei']
            elif par == 'RawIntDen':
                nucleus_component = derived_values['RawIntDen_nuclei']
            else:
                nucleus_component = np.nan
                
            if par == 'Area':
                cytoplasm_component = derived_values['area_cytoplasm']
            elif par == 'Mean':
                cytoplasm_component = derived_values['MGV_cytoplasm']
            elif par == 'StdDev':
                cytoplasm_component = derived_values['StdDev_cytoplasm']
            elif par == 'IntDen':
                cytoplasm_component = derived_values['ID_cytoplasm']
            else:
                cytoplasm_component = np.nan
                
            nuc_cell_ratio_name = f'{mapped_name}_nuclei_cell_ratio'
            derived_values[nuc_cell_ratio_name] = safe_ratio(nucleus_component, cell_component)
            
            if par != 'RawIntDen':
                nuc_cyto_ratio_name = f'{mapped_name}_nuclei_cytoplasm_ratio'
                derived_values[nuc_cyto_ratio_name] = safe_ratio(nucleus_component, cytoplasm_component)

        for par in parameters:
            all_data_by_condition[condition][par].append(cell_data.get(par, np.nan))
        
        for par in all_derived_metrics:
            all_data_by_condition[condition][par].append(derived_values[par])

    field_cell_counts[folder] = processed_cells_in_field


print("\nConsolidating and saving data...")
all_parameters_to_save = parameters + all_derived_metrics

for par in all_parameters_to_save:
    if not par.strip():
        continue
        
    data_for_df = {
        cond: pd.Series(all_data_by_condition[cond][par])
        for cond in all_conditions
    }
    df_par = pd.DataFrame(data_for_df)

    df_par = df_par.dropna(how='all')

    if par in PARAMETER_MAP:
        output_name = PARAMETER_MAP[par]
    elif par.startswith('#'):
        output_name = par.replace("#", "Number_")
    else:
        output_name = par
        
    csv_file = os.path.join(output_dir, f'combined_{output_name}.csv')
    df_par.to_csv(csv_file, index=False)
    print(f"Saved {par} data to: {csv_file}")


print("\nSaving all individual nucleus RAW data (per nucleus, by parameter)...")

for param_raw, data_by_cond in single_nuclei_data_by_parameter.items():
    if not any(data_by_cond.values()):
        continue
        
    data_by_cond_series = {
        cond: pd.Series(data_by_cond[cond])
        for cond in all_conditions
    }
    df_par_single = pd.DataFrame(data_by_cond_series)

    mapped_name = PARAMETER_MAP.get(param_raw, param_raw)
    output_name = f"{mapped_name}_single_nuclei"
    
    csv_file = os.path.join(output_dir, f'combined_{output_name}.csv')
    df_par_single.to_csv(csv_file, index=False)
    print(f"Saved Single Nucleus {param_raw} data to: {csv_file}")

print("\nSaving individual nucleus Eop data (per nucleus)...")
if any(single_Eop_data_by_condition.values()):
    data_by_cond_series = {
        cond: pd.Series(single_Eop_data_by_condition[cond])
        for cond in all_conditions
    }
    df_eop_single = pd.DataFrame(data_by_cond_series)
    
    output_name = "Eop_single_nuclei"
    csv_file = os.path.join(output_dir, f'combined_{output_name}.csv')
    df_eop_single.to_csv(csv_file, index=False)
    print(f"Saved Single Nucleus Eop data to: {csv_file}")
else:
    print("Skipped Eop calculation for Single Nuclei: No nucleus data found.")


print("\nCalculating and saving CoV data...")

data_mean_cell_series = {cond: pd.Series(all_data_by_condition[cond]['Mean']) for cond in all_conditions}
data_stddev_cell_series = {cond: pd.Series(all_data_by_condition[cond]['StdDev']) for cond in all_conditions}

df_mean_cell = pd.DataFrame(data_mean_cell_series)
df_stddev_cell = pd.DataFrame(data_stddev_cell_series)

df_cov_cell = pd.DataFrame()

for cond in all_conditions:
    cov_series = df_stddev_cell[cond] / df_mean_cell[cond]
    df_cov_cell[cond] = cov_series.replace([np.inf, -np.inf], np.nan)

df_cov_cell = df_cov_cell.dropna(how='all')
csv_file_cell_cov = os.path.join(output_dir, 'combined_CoV_cell.csv')
df_cov_cell.to_csv(csv_file_cell_cov, index=False)
print(f"Saved CoV for Cell to: {csv_file_cell_cov}")


data_mean_cyto_series = {cond: pd.Series(all_data_by_condition[cond]['MGV_cytoplasm']) for cond in all_conditions}
data_stddev_cyto_series = {cond: pd.Series(all_data_by_condition[cond]['StdDev_cytoplasm']) for cond in all_conditions}

df_mean_cyto = pd.DataFrame(data_mean_cyto_series)
df_stddev_cyto = pd.DataFrame(data_stddev_cyto_series)

df_cov_cyto = pd.DataFrame()

for cond in all_conditions:
    cov_series = df_stddev_cyto[cond] / df_mean_cyto[cond]
    df_cov_cyto[cond] = cov_series.replace([np.inf, -np.inf], np.nan)

df_cov_cyto = df_cov_cyto.dropna(how='all')
csv_file_cyto_cov = os.path.join(output_dir, 'combined_CoV_cytoplasm.csv')
df_cov_cyto.to_csv(csv_file_cyto_cov, index=False)
print(f"Saved CoV for Cytoplasm to: {csv_file_cyto_cov}")


cov_nuclei_results = {}
for cond in all_conditions:
    mean_data = pd.Series(single_nuclei_data_by_parameter['Mean'][cond])
    stddev_data = pd.Series(single_nuclei_data_by_parameter['StdDev'][cond])
    
    if not mean_data.empty:
        cov_series = stddev_data / mean_data
        cov_nuclei_results[cond] = cov_series.replace([np.inf, -np.inf], np.nan)

if cov_nuclei_results:
    df_cov_nuclei = pd.DataFrame(cov_nuclei_results)
    df_cov_nuclei = df_cov_nuclei.dropna(how='all')
    csv_file_nuclei_cov = os.path.join(output_dir, 'combined_CoV_nuclei_single.csv')
    df_cov_nuclei.to_csv(csv_file_nuclei_cov, index=False)
    print(f"Saved CoV for Single Nuclei to: {csv_file_nuclei_cov}")
else:
    print("Skipped CoV calculation for Single Nuclei: No nucleus data found.")

print("\nSaving Cells per Field data in quantifiable (wide) format...")

counts_by_condition = {cond: [] for cond in all_conditions}
fields_processed = sorted(field_cell_counts.keys())

for folder in fields_processed:
    condition = folder_map[folder]
    count = field_cell_counts[folder]
    counts_by_condition[condition].append(count)

data_for_df = {
    cond: pd.Series(counts_by_condition[cond])
    for cond in all_conditions
}

df_field_counts = pd.DataFrame(data_for_df)
df_field_counts = df_field_counts.dropna(how='all')

csv_file_field_counts = os.path.join(output_dir, 'cells_per_field.csv')
df_field_counts.to_csv(csv_file_field_counts, index=False)
print(f"Saved processed cell counts per field (wide format) to: {csv_file_field_counts}")


print("\nData consolidation complete. The combined data is in the 'combined_output' directory.")