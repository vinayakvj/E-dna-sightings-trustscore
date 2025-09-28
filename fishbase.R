suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(arrow)
})

in_path <- "C://Users//rosam//Downloads//all_voyages.csv//all_voyages.csv" #ADD YOUR FILE

cache_dir <- "fishbase_cache"
dir.create(cache_dir, showWarnings = FALSE)

dl <- function(url, fname) {
  path <- file.path(cache_dir, fname)
  if (!file.exists(path)) download.file(url, path, mode = "wb", quiet = TRUE)
  path
}

species_path <- dl(
  "https://huggingface.co/datasets/cboettig/fishbase/resolve/main/data/fb/v24.07/parquet/species.parquet?download=true",
  "species.parquet"
)
syn_path <- dl(
  "https://huggingface.co/datasets/cboettig/fishbase/resolve/main/data/fb/v24.07/parquet/synonyms.parquet?download=true",
  "synonyms.parquet"
)
country_path <- dl(
  "https://huggingface.co/datasets/cboettig/fishbase/resolve/main/data/fb/v24.07/parquet/country.parquet?download=true",
  "country.parquet"
)

species_df <- read_parquet(species_path) %>% select(SpecCode, Genus, Species)
syn_df     <- read_parquet(syn_path)     %>% select(SpecCode, SynGenus, SynSpecies)

country_raw <- read_parquet(country_path)
spec_col   <- grep("spec", names(country_raw), ignore.case = TRUE, value = TRUE)[1]
ccode_col  <- grep("c[_]?code", names(country_raw), ignore.case = TRUE, value = TRUE)[1]

country_df <- country_raw %>%
  select(SpecCode = all_of(spec_col),
         C_Code   = all_of(ccode_col))

au_df <- country_df %>%
  filter(C_Code %in% c(36,"36","036"))

is_in_au <- function(code) {
  if (is.na(code)) return(NA)
  any(au_df$SpecCode == code)
}

`%||%` <- function(a,b) if (is.null(a) || is.na(a) || (is.character(a) && trimws(a)=="")) b else a

name_to_speccode <- function(sciname) {
  sciname <- str_squish(sciname %||% "")
  if (sciname == "" || !str_detect(sciname, "\\s")) return(NA_integer_)
  parts <- str_split(sciname, "\\s+", n = 2, simplify = TRUE)
  G <- parts[1,1]; S <- parts[1,2]
  hit <- species_df %>% filter(Genus == G, Species == S) %>% slice_head(n = 1)
  if (nrow(hit) == 1) return(hit$SpecCode[[1]])
  syn <- syn_df %>% filter(SynGenus == G, SynSpecies == S) %>% slice_head(n = 1)
  if (nrow(syn) == 1) return(syn$SpecCode[[1]])
  NA_integer_
}

speccode_to_name <- function(code) {
  if (is.na(code)) return(NA_character_)
  r <- species_df %>% filter(SpecCode == code) %>% slice_head(n = 1)
  if (nrow(r) == 0) return(NA_character_)
  paste(r$Genus, r$Species)
}

df <- read_csv(in_path, show_col_types = FALSE)

process_row <- function(row) {
  row <- as.list(row)
  sp  <- as.character(row$species %||% NA)
  lca <- as.character(row$Species_In_LCA %||% NA)
  asv <- as.character(row$ASV %||% NA)
  samp <- as.character(row$sample %||% NA)
  
  if (!is.na(sp) && tolower(str_squish(sp)) != "dropped") {
    sc <- name_to_speccode(sp)
    res <- tibble::tibble(
      ASV            = asv,
      sample         = samp,
      case           = "exact_species",
      input_name     = sp,
      candidate_name = speccode_to_name(sc),
      SpecCode       = sc,
      in_AU          = is_in_au(sc)
    )
    res <- res %>% mutate(AU_summary = ifelse(is.na(in_AU), NA, ifelse(in_AU, "T", "F")))
    return(res)
  } else if (!is.na(lca) && lca != "") {
    cand_list <- str_split(lca, ",")[[1]] %>% str_squish()
    scs <- map_int(cand_list, name_to_speccode)
    vals <- map_lgl(scs, is_in_au)
    summary_str <- paste(ifelse(is.na(vals), NA, ifelse(vals, "T", "F")), collapse = ",")
    tibble::tibble(
      ASV            = asv,
      sample         = samp,
      case           = "lca_candidates_list",
      input_name     = lca,
      candidate_name = map_chr(scs, speccode_to_name),
      SpecCode       = scs,
      in_AU          = vals,
      AU_summary     = summary_str
    )
  } else {
    tibble::tibble(
      ASV            = asv,
      sample         = samp,
      case           = "no_species_info",
      input_name     = NA_character_,
      candidate_name = NA_character_,
      SpecCode       = NA_integer_,
      in_AU          = NA,
      AU_summary     = NA_character_
    )
  }
}

out <- df %>%
  mutate(.row_id = row_number()) %>%
  split(.$.row_id) %>%
  map_dfr(process_row)

cat("\n--- Candidate species per row (with Australia flags) ---\n")
print(out, n = 40)

cat("\n--- Summary per input row (all candidates concatenated) ---\n")
summary_per_input <- out %>%
  group_by(ASV, sample, case, input_name, AU_summary) %>%
  summarise(n_candidates = n(), .groups = "drop")
print(summary_per_input, n = 40)



