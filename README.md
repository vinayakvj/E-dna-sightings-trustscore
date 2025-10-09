# Species Evidence & Distance Pipeline (R)

## Run
1. Open `main.R` and set `file_path` to your local CSV path.
2. Install packages:
   ```r
   pkgs <- c("dplyr","readr","stringr","tidyr","tibble","purrr","rlang",
             "sf","geosphere","robis","worrms","arrow")
   install.packages(setdiff(pkgs, rownames(installed.packages())))
