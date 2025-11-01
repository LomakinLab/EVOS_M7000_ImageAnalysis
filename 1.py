import os
import shutil
from pathlib import Path

def organize_channels(root_dir: Path):
    dapi_folder_name = "DAPI_Channel"
    actin_folder_name = "Actin_Channel"
    poi1_folder_name = "PoI1_Channel"
    poi2_folder_name = "PoI2_Channel"
    
    dapi_folder = root_dir / dapi_folder_name
    actin_folder = root_dir / actin_folder_name
    poi1_folder = root_dir / poi1_folder_name
    poi2_folder = root_dir / poi2_folder_name
    
    dapi_folder.mkdir(exist_ok=True)
    actin_folder.mkdir(exist_ok=True)
    poi1_folder.mkdir(exist_ok=True)
    poi2_folder.mkdir(exist_ok=True)
    
    print(f"--- Channel Organization Started ---")
    print(f"Output folders created:")
    print(f"- {dapi_folder.name}")
    print(f"- {actin_folder.name}")
    print(f"- {poi1_folder.name}")
    print(f"- {poi2_folder.name}")

    DAPI_PATTERN = "*_R_p00_0_A*d0.TIF"
    ACTIN_PATTERN = "*_R_p00_0_A*d2.TIF"
    POI1_PATTERN = "*_R_p00_0_A*d1.TIF"
    POI2_PATTERN = "*_R_p00_0_A*d3.TIF"
    
    all_channel_names = (dapi_folder_name, actin_folder_name, poi1_folder_name, poi2_folder_name)

    for folder in root_dir.iterdir():
        
        if folder.is_dir() and folder.name not in all_channel_names:
            
            print(f"\nScanning source folder: {folder.name}/")
            
            dapi_files = list(folder.glob(DAPI_PATTERN))
            for file in dapi_files:
                shutil.copy2(file, dapi_folder / file.name)
                print(f"  -> Copied DAPI file (d0): {file.name}")
            
            if not dapi_files:
                print("  -> No DAPI (d0.TIF) files found matching the pattern.")

            actin_files = list(folder.glob(ACTIN_PATTERN))
            for file in actin_files:
                shutil.copy2(file, actin_folder / file.name)
                print(f"  -> Copied Actin file (d2): {file.name}")
                
            if not actin_files:
                print("  -> No Actin (d2.TIF) files found matching the pattern.")

            poi1_files = list(folder.glob(POI1_PATTERN))
            for file in poi1_files:
                shutil.copy2(file, poi1_folder / file.name)
                print(f"  -> Copied PoI1 file (d1): {file.name}")
                
            if not poi1_files:
                print("  -> No PoI1 (d1.TIF) files found matching the pattern.")

            poi2_files = list(folder.glob(POI2_PATTERN))
            for file in poi2_files:
                shutil.copy2(file, poi2_folder / file.name)
                print(f"  -> Copied PoI2 file (d3): {file.name}")
                
            if not poi2_files:
                print("  -> No PoI2 (d3.TIF) files found matching the pattern.")

    print("\n--- Organization Complete! ---")
    print(f"Total DAPI files copied: {len(list(dapi_folder.iterdir()))}")
    print(f"Total Actin files copied: {len(list(actin_folder.iterdir()))}")
    print(f"Total PoI1 files copied: {len(list(poi1_folder.iterdir()))}")
    print(f"Total PoI2 files copied: {len(list(poi2_folder.iterdir()))}")

if __name__ == "__main__":
    root_directory = Path.cwd()
    organize_channels(root_directory)