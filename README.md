# Species Evidence & Distance Pipeline (R)

## Run
1. Open `main.R` and set `file_path` to your local CSV path.
2. Install packages:
   ```r
   pkgs <- c("dplyr","readr","stringr","tidyr","tibble","purrr","rlang",
             "sf","geosphere","robis","worrms","arrow")
   install.packages(setdiff(pkgs, rownames(installed.packages())))



## 📁 Data Folder (Pre-Zipped Outputs)

All project results and cached OBIS data are stored inside the `data/` folder.  
This ensures the entire pipeline can be reviewed or rerun instantly without waiting for OBIS downloads.

### 1️⃣ species_in_LCA_truth_per_row.zip
This ZIP contains all four core output CSVs:
- `data_OBIS_enriched_distances.csv`  
- `OBIS_distance_quantiles_per_species.csv`  
- `species_in_LCA_truth_table.csv`  
- `species_in_LCA_truth_per_row.csv`

Each of these files contains the processed outputs of the OBIS validation pipeline:
- **data_OBIS_enriched_distances.csv:** Row-level proximity and spatial flags  
- **OBIS_distance_quantiles_per_species.csv:** Per-species distance statistics  
- **species_in_LCA_truth_table.csv:** Species-level AU presence evidence  
- **species_in_LCA_truth_per_row.csv:** Row × species evidence mapping  

###  SpeciesPoints_OBIS.zip
This ZIP contains the full cached OBIS species occurrence data.  
It corresponds to the `Data/SpeciesPoints_OBIS/` folder that the script generates automatically.  
Having it zipped here means others can run `main.R` without re-downloading OBIS data.

### How to Use
1. Go to the `data/` folder on GitHub.  
2. Click on each ZIP file → **Download**.  
3. Unzip both files in your project directory.  
   After unzipping, your structure should look like this:
